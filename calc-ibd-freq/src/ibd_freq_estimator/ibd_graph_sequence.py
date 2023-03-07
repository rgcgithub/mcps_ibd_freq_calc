import heapq
import numpy as np
import scipy.sparse


class IBDGraphSequence:
    """Class to interate through a sequence of graphs defined by IBD segments.

    IBD graphs along the genome are defined by a sequencing of breakpoints that
    determine when the structure if the IBD graph changes. Consider the
    following example of IBD segments along the genome:

        IBD SEG         ========   ========
        IBD SEG     =======        =====      ===  ====
         GENOME  ***************************************
    BREAKPOINTS     ^   ^  ^    ^  ^    ^  ^  ^  ^ ^   ^

    At each breakpoint, depicted by a ^, the structure of the IBD graph changes.
    We want to (1) find each breakpoint, then (2) interrogate the graph in the
    span between two breakpoints.

    For each iteration, we maintain a variable defining the genomic span ---
    the start and end points of the current IBD graph --- and a variable
    with the collection of current IBD segments that overlap that span.

    (1) Finding breakpoints
        When we observe a new IBD segment, there are two cases:
            A) The start position is the same as the current breakpoint. In
               this case we add the IBD segment to the current graph. We
               recompute the next breakpoint to be the smallest of the
               endpoints of the current IBD segments.

            B) The start position is greater than the current breakpoint. We
               set the next breakpoint to the smaller of the start position
               of the IBD segment and the stored value of the next breakpoint.

        We add IBD segments until we hit case B. This gives us the genomic span
        of the current IBD graph.

    (2) Advancing to the next breakpoint
        Given the current set of IBD segments, we advance the current breakpoint
        to the next breakpoint, and set the next breakpoint to the smallest
        of the current IBD segments. We drop IBD segments that no longer overlap
        with the current span. Then we repeat (1) to get the end point of
        the current genomic span.

    We iterate through (1) and (2) until we have advanced across all graphs
    defined by the IBD segments.
    """

    def __init__(self):
        # matrix of start_pos, end_pos, sample-hap ID 1, sample-hap ID 2
        self.active_ibd_segs = np.zeros((10000, 4))
        self.incoming_ibd_segs = np.zeros(self.active_ibd_segs.shape)
        self.active_ibd_segs_last_idx = 0
        self.incoming_ibd_segs_last_idx = 0

        self.current_breakpoint = 0
        self.next_breakpoint = 1

        # sparse IBD graph
        self.data = [] # list of nonzero edge weights
        self.idxes = [] # list of tuples of nonzero edges in IBD graph

        # map sample ids + haplotype ids to integer indexes
        self.sid_hid_to_idx = {}
        # integer indexes to sample id + haplotype id
        self.idx_to_sid_hid = {}
        # integer indexes of active samples in the IBD graph
        self.active_idx = set()


    def add_ibd_segment(self, ibd_segment):
        """Checks if an IBD segment overlaps the current IBD graph, and if so
        updates the graph. Otherwise the IBD segment is added to the list of
        incoming IBD segments.
        """
        if self.overlaps_current_ibd_graph(ibd_segment):
            self.add_active_ibd_segment(ibd_segment)
            updates_ibd_graph = True
            return updates_ibd_graph
        else:
            self.add_incoming_ibd_segment(ibd_segment)
            updates_ibd_graph = False
            return updates_ibd_graph


    def add_active_ibd_segment(self, ibd_segment):
        """Adds the IBD segment to the current IBD graph. Raises an error
        if the start position of the IBD segment is greater than the current
        breakpoint.
        """
        assert ibd_segment.start_pos == self.current_breakpoint
        self.next_breakpoint = min(self.next_breakpoint, ibd_segment.end_pos + 1)

        # update IBD graph
        idx1 = self._get_sample_hap_idx(
            ibd_segment.sid1,
            ibd_segment.hid1
        )
        idx2 = self._get_sample_hap_idx(
            ibd_segment.sid2,
            ibd_segment.hid2
        )

        self.active_ibd_segs[self.active_ibd_segs_last_idx] = [
            ibd_segment.start_pos,
            ibd_segment.end_pos,
            idx1,
            idx2
        ]
        self.active_ibd_segs_last_idx += 1
        if self.active_ibd_segs_last_idx == self.active_ibd_segs.shape[0]:
            self.active_ibd_segs = self._resize_arr(self.active_ibd_segs)

        self.active_idx.add(idx1)
        self.active_idx.add(idx2)

        self.data.append(1)
        self.idxes.append([idx1, idx2])


    def add_incoming_ibd_segment(self, ibd_segment):
        """Check the next breakpoint against the IBD segment. Reduces the
        next breakpoint to min(next_breakpoint, ibd_segment.start_pos). Adds
        the IBD segment to the list of IBD segments after the current graph.
        """
        m = min(ibd_segment.start_pos, self.next_breakpoint)
        self.next_breakpoint = m

        idx1 = self._get_sample_hap_idx(
            ibd_segment.sid1,
            ibd_segment.hid1
        )
        idx2 = self._get_sample_hap_idx(
            ibd_segment.sid2,
            ibd_segment.hid2
        )

        self.incoming_ibd_segs[self.incoming_ibd_segs_last_idx] = [
            ibd_segment.start_pos,
            ibd_segment.end_pos,
            idx1,
            idx2
        ]
        self.incoming_ibd_segs_last_idx += 1

        if self.incoming_ibd_segs_last_idx == self.incoming_ibd_segs.shape[0]:
            self.incoming_ibd_segs = self._resize_arr(self.incoming_ibd_segs)


    def overlaps_current_ibd_graph(self, ibd_segment):
        """Returns True if the IBD segments overlaps the span of the
        current IBD graph.
        """
        return ibd_segment.start_pos == self.current_breakpoint


    def advance_to_next_breakpoint(self):
        """Advances the graph to the next breakpoint. Drops IBD segments that
        no longer overlap the genomic span of the IBD graph.
        """
        current_size = self.active_ibd_segs.shape[0]

        self.current_breakpoint = self.next_breakpoint
        self.next_breakpoint = np.inf

        keep_idx = \
            self.active_ibd_segs[:self.active_ibd_segs_last_idx, 1] \
            >= self.current_breakpoint
        n_keep = keep_idx.sum()
        self.active_ibd_segs[:n_keep] = \
            self.active_ibd_segs[:self.active_ibd_segs_last_idx][keep_idx]
        self.active_ibd_segs_last_idx = n_keep
        assert self.active_ibd_segs_last_idx == 0 \
            or np.all(self.active_ibd_segs[:self.active_ibd_segs_last_idx, 0] <= self.current_breakpoint)

        new_active_idx = np.searchsorted(
            self.incoming_ibd_segs[:self.incoming_ibd_segs_last_idx, 0],
            self.current_breakpoint,
            side="right"
        )
        new_active_ibd_segs = self.incoming_ibd_segs[:new_active_idx]
        self.incoming_ibd_segs = self.incoming_ibd_segs[new_active_idx:]
        self.incoming_ibd_segs_last_idx -= new_active_idx
        assert np.all(new_active_ibd_segs[:,0] == self.current_breakpoint)
        assert self.incoming_ibd_segs_last_idx == 0 or \
            np.all(self.incoming_ibd_segs[:self.incoming_ibd_segs_last_idx, 0] > self.current_breakpoint)

        add_rows = max(
            0,
            current_size - self.active_ibd_segs_last_idx - new_active_ibd_segs.shape[0]
        )
        stack = []
        if self.active_ibd_segs_last_idx > 0:
            stack.append(self.active_ibd_segs[:self.active_ibd_segs_last_idx])
        if new_active_ibd_segs.shape[0] > 0:
            stack.append(new_active_ibd_segs)
        # combining segments would not trigger a resize
        if add_rows > 0:
            stack.append(np.zeros((add_rows, 4)))
        else:
            stack.append(np.zeros((current_size, 4)))
        self.active_ibd_segs = np.vstack(stack)

        self.active_ibd_segs_last_idx = self.active_ibd_segs_last_idx + new_active_ibd_segs.shape[0]

        potential_breakpoint = np.inf
        if self.active_ibd_segs_last_idx > 0:
            potential_breakpoint = self.active_ibd_segs[
                :self.active_ibd_segs_last_idx,
                1
            ].min() + 1
        if self.incoming_ibd_segs_last_idx > 0:
            potential_breakpoint = min(
                potential_breakpoint,
                self.incoming_ibd_segs[
                    :self.incoming_ibd_segs_last_idx,
                    0
                ].min()
            )
        self.next_breakpoint = min(
            self.next_breakpoint,
            potential_breakpoint
        )

        assert self.current_breakpoint < self.next_breakpoint or \
            (np.isinf(self.current_breakpoint) and np.isinf(self.next_breakpoint))

        self.active_idx = set(
            self.active_ibd_segs[:self.active_ibd_segs_last_idx, 2]
        ).union(
            self.active_ibd_segs[:self.active_ibd_segs_last_idx, 3]
        )

        # this operation is expensive, so let's not do this every time
        if len(self.active_idx) > 5*len(self.sid_hid_to_idx):
            self._reindex_samples()
        else:
            self.idxes = list(
                self.active_ibd_segs[:self.active_ibd_segs_last_idx,2:]
            )
            self.data = [1 for i in range(self.active_ibd_segs_last_idx)]


    def _reindex_samples(self):
        """Reset the integer index for sample-haplotype pairs and
        drop sample-haplotypes that are no longer part of the incoming
        or active IBD segments.
        """
        active_sample_haps = []
        for row in self.active_ibd_segs:
            active_sample_haps.append(
                (
                    self.get_sample_hap[row[2]],
                    self.get_sample_hap[row[3]]
                )
            )

        incoming_sample_haps = []
        for row in self.incoming_ibd_segs:
            incoming_sample_haps.append(
                (
                    self.get_sample_hap[row[2]],
                    self.get_sample_hap[row[3]]
                )
            )

        # update IBD graph
        # map sample ids + haplotype ids to integers
        self.sid_hid_to_idx = {}
        self.idx_to_sid_hid = {}
        self.active_idx = set()

        sid_hid1 = []
        sid_hid2 = []
        data = []
        idxes = []
        for row in active_sample_haps:
            sid_hid1.append(
                self._get_sample_hap_idx(row[0])
            )
            sid_hid2.append(
                self._get_sample_hap_idx(row[1])
            )
            idxes.append(
                [sid_hid1[-1], sid_hid2[-1]]
            )
            data.append(1)
        self.active_ibd_segs[:,2] = sid_hid1
        self.active_ibd_segs[:,3] = sid_hid2
        self.data = data
        self.idxes = idxes

        sid_hid1 = []
        sid_hid2 = []
        for row in incoming_sample_haps:
            sid_hid1.append(
                self._get_sample_hap_idx(row[0])
            )
            sid_hid2.append(
                self._get_sample_hap_idx(row[1])
            )
        self.incoming_ibd_segs[:,2] = sid_hid1
        self.incoming_ibd_segs[:,3] = sid_hid2


    def get_connected_components(self):
        """Return connected components of the IBD graph with size > 1.

        Returns
        -------
            conn_comps : list[set]
                A list of connected components represented by the set
                of integer ids of that component.

        Notes
        -----
            Transform sample ids + haplotype ids to integer ids with
            get_sample_hap_idx(sid, hid)
        """
        idxes = np.array(self.idxes)
        N = int(
            max(self.active_idx) + 1 if len(self.active_idx) > 0 else 0
        )

        if N == 0:
            return []

        csr = scipy.sparse.csr_matrix(
            (
                self.data,
                (idxes[:,0], idxes[:,1]),
            ),
            shape=(N, N)
        )

        n_components, labels = scipy.sparse.csgraph.connected_components(
            csr,
            directed=False,
            return_labels=True
        )

        label_to_tuple = {}
        for idx,l in enumerate(labels):
            if not l in label_to_tuple:
                label_to_tuple[l] = set()

            label_to_tuple[l].add(idx)

        conn_comps = [
            label_to_tuple[l] for l in label_to_tuple \
            if len(label_to_tuple[l]) > 1
        ]
        return conn_comps


    def get_sample_hap_idx(self, sid, hid):
        """Returns the integer id of the sample + haplotype id if it is part of
        the current IBD graph, otherwise returns -1
        """
        idx =  self.sid_hid_to_idx.get(sid + "." + hid, -1)
        if idx == -1 or not idx in self.active_idx:
            return -1
        else:
            return idx


    def has_at_least_one_edge(self):
        return self.active_ibd_segs_last_idx > 0


    def has_incoming_ibd_segments(self):
        return self.incoming_ibd_segs_last_idx > 0


    def has_active_ibd_segments(self):
        return self.has_at_least_one_edge()


    def get_span(self):
        current_breakpoint = int(self.current_breakpoint) \
            if np.isfinite(self.current_breakpoint) \
            else np.inf
        next_breakpoint = int(self.next_breakpoint) \
            if np.isfinite(self.next_breakpoint) \
            else np.inf

        return current_breakpoint, next_breakpoint - 1


    def _get_sample_hap_idx(self, sid, hid):
        """Transform a sample id + haplotype id to an integer id.
        """
        sid_hid = sid + "." + hid
        if not sid_hid in self.sid_hid_to_idx:
            self.sid_hid_to_idx[sid_hid] = len(self.sid_hid_to_idx)

        idx = self.sid_hid_to_idx[sid_hid]
        self.idx_to_sid_hid[idx] = sid_hid

        return idx


    def _resize_arr(self, arr):
        N = arr.shape[0]
        zeros = np.zeros((N, 4))
        return np.vstack((arr, zeros))


    def get_sample_hap(self, idx):
        """Transform a integer id of a sample + haplotype to the tuple
        (sample_id, haplotype).
        """
        sid_hid = self.idx_to_sid_hid[idx]
        return tuple(sid_hid.split("."))
