import errno
import numpy as np
import scipy
import logging
import os
import sys
import time
import xopen

from ibd_freq_estimator.ibd_graph_sequence import IBDGraphSequence

import ibd_freq_estimator.logging_config
logger = logging.getLogger(__name__)


class IBDSegment:

    def __init__(self, sid1, sid2, hid1, hid2, chrom, start_pos, end_pos):
        self.sid1 = sid1
        self.sid2 = sid2
        self.hid1 = hid1
        self.hid2 = hid2
        self.chrom = chrom
        self.start_pos = start_pos
        self.end_pos = end_pos


    def __str__(self):
        return f"IBDSegment(sid1: {self.sid1}, sid2: {self.sid2}, hid1: {self.hid1}, hid2: {self.hid2}, chr: {self.chrom}, start: {self.start_pos}, end: {self.end_pos})"


    def __repr__(self):
        return self.__str__()


class IBDGraph:
    """Wrapper around IBDGraphSequence. Prevent exposing internal methods to
    users of the IBDFile class.
    """

    def __init__(self, ibd_graph_sequence):
        self.ibd_graph_sequence = ibd_graph_sequence


    def get_span(self):
        """Returns the start and end coordinates (inclusive) along the genome
        of the current IBD graph.
        """
        return self.ibd_graph_sequence.get_span()


    def get_connected_components(self):
        """Find the connected components of the current IBD graph.

        Returns
        -------
            conn_comps : list[set[int]]
                A collection of sets of connected components with samples
                specified by their integer id.
        """
        conn_comp = self.ibd_graph_sequence.get_connected_components()
        return conn_comp


    def get_sample_hap_idx(self, sid, hid):
        """Returns the integer id of hap 'hap' in sample 'sid'. Returns
        -1 if the sample is not part of the current IBD graph.
        """
        return self.ibd_graph_sequence.get_sample_hap_idx(sid, hid)


    def get_sample_hap(self, idx):
        """Returns (sample_id, haplotype_id) from a integer id in the
        IBD graph.
        """
        return self.ibd_graph_sequence.get_sample_hap(idx)


    def has_at_least_one_edge(self):
        return self.ibd_graph_sequence.has_at_least_one_edge()


    def size(self):
        return self.ibd_graph_sequence.active_ibd_segs_last_idx


class InvalidIBDFileError(Exception):
    pass


class IBDFile:
    """Iterate through IBD graphs defined by and IBD calls. IBD calls must be
    sorted by start position.

    Parameters
    ----------
        path : str
            Path to IBD file
        sid1_col : int
            Zero-based column index of the first sample id
        sid2_col : int
            Zero-based column index of the second sample id
        hid1_col : int
            Zero-based column index of the haplotype id of the first sample
        hid2_col : int
            Zero-based column index of the heplotype id of the second sample
        chr_col : int
            Zero-based column index of the chromosome of the IBD segment
        start_pos_col : int
            Zero-based column index of the starting position of the IBD segment
        end_pos_col : int
            Zero-based column index of the ending position of the IBD segment
        has_header : bool
            True if the IBD file has a header line
        hap1_symbol : str
            Label for haplotype 1 in the IBD file (e.g. 0, 1)
        hap2_symbol : str
            Label for haplotype 2 in the IBD file (e.g. 1, 2)
    """

    def __init__(
        self,
        path,
        sid1_col,
        sid2_col,
        hid1_col,
        hid2_col,
        chr_col,
        start_pos_col,
        end_pos_col,
        has_header,
        hap1_symbol = None,
        hap2_symbol = None
    ):
        if not os.path.exists(path):
            raise FileNotFoundError(
                eerno.ENOENT,
                os.strerror(errno.ENOENT),
                path
            )

        self.path = path
        self.sid1_col = sid1_col
        self.sid2_col = sid2_col
        self.hid1_col = hid1_col
        self.hid2_col = hid2_col
        self.chr_col = chr_col
        self.start_pos_col = start_pos_col
        self.end_pos_col = end_pos_col
        self.handle = xopen.xopen(self.path)
        self.ibd_graph = IBDGraphSequence()
        self.eof = False
        self.last_start_pos = -1

        # symbols denoting haplotype 1 and haplotype 2
        if hap1_symbol is None or hap2_symbol is None:
            self.hap1_symbol, self.hap2_symbol = self.get_haplotype_symbols(
                path,
                hid1_col,
                hid2_col,
                has_header
            )
        else:
            self.hap1_symbol = hap1_symbol
            self.hap2_symbol = hap2_symbol

        if has_header:
            self.handle.readline()

        ibd_seg = self._read_ibd_segment()
        if ibd_seg != "":
            self.ibd_graph.add_ibd_segment(ibd_seg)


    def __iter__(self):
        return self


    def __next__(self):
        if self.eof and \
            not self.ibd_graph.has_active_ibd_segments() and \
            not self.ibd_graph.has_incoming_ibd_segments():
            self.handle.close()
            raise StopIteration

        self.ibd_graph.advance_to_next_breakpoint()
        updates_ibd_graph = True
        while updates_ibd_graph:
            ibd_seg = self._read_ibd_segment()

            if ibd_seg == "":
                updates_ibd_graph = False
                self.eof = True
                break

            if ibd_seg.start_pos < self.last_start_pos:
                raise ValueError("IBD segments must be sorted by start position")
            else:
                self.last_start_pos = ibd_seg.start_pos
                updates_ibd_graph = self.ibd_graph.add_ibd_segment(ibd_seg)

        ibd_graph = IBDGraph(self.ibd_graph)
        return ibd_graph


    def __exit__(self, *args):
        self.handle.close()


    def _read_ibd_segment(self):
        line = self.handle.readline()
        if line != "":
            line = line.strip("\n").split()
            ibd_seg = IBDSegment(
                line[self.sid1_col],
                line[self.sid2_col],
                line[self.hid1_col],
                line[self.hid2_col],
                line[self.chr_col],
                int(line[self.start_pos_col]),
                int(line[self.end_pos_col])
            )
            return ibd_seg
        else:
            return ""


    def close(self):
        self.handle.close()


    def get_haplotype_symbols(self, path, hid1_col, hid2_col, has_header):
        """Find the symbols identifying haplotype 1 and haplotype 2.
        """
        hsymbols = set()
        lines_read = 0
        with xopen.xopen(path) as f:
            for i,line in enumerate(f):
                if has_header and i == 0: continue
                line = line.strip("\n").split()
                hsymbols.add(line[hid1_col])
                hsymbols.add(line[hid2_col])
                lines_read += 1

                # scan the file until we find all haplotype symbols
                if i % 100000 == 0 and len(hsymbols) >= 2:
                    break

        if len(hsymbols) != 2 and lines_read > 0:
            raise InvalidIBDFileError(f"found {len(hsymbols)} haplotype symbols ({hsymbols}): expected 2")
        elif lines_read == 0:
            logger.warning("ibd file is empty")
            hsymbols.add("1")
            hsymbols.add("2")

        hsymbols = sorted(list(hsymbols))
        return hsymbols[0], hsymbols[1]
