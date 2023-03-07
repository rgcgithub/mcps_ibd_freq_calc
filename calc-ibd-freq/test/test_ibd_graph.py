import unittest

from src.ibd_freq_estimator.ibd_graph_sequence import IBDGraphSequence
from src.ibd_freq_estimator.ibd_file import IBDSegment
from test.test_graphs import get_test_graph1

class TestIBDGraphSequence(unittest.TestCase):

    def test_breakpoints(self):
        test_ibd_file = "./test/data/ibd_segs.txt"
        ibd_segs = load_ibd_seg_file(test_ibd_file)

        ibd_graph = IBDGraphSequence()
        breakpoints = []

        conn_comps, spans = get_test_graph1()
        true_breakpoints = [s[0] for s in spans]

        for ibd_seg in ibd_segs:
            updates_ibd_graph = ibd_graph.add_ibd_segment(ibd_seg)
            if not updates_ibd_graph:
                ibd_graph.advance_to_next_breakpoint()
                breakpoints.append(ibd_graph.current_breakpoint)

        while ibd_graph.has_incoming_ibd_segments() or ibd_graph.has_active_ibd_segments():
            ibd_graph.advance_to_next_breakpoint()
            breakpoints.append(ibd_graph.current_breakpoint)

        for b in true_breakpoints:
            self.assertTrue(b in breakpoints, msg=f"missing {b} from breakpoints")


    def test_graph_construction(self):
        test_ibd_file = "./test/data/ibd_segs.txt"
        ibd_segs = load_ibd_seg_file(test_ibd_file)

        true_conn_comps, true_spans = get_test_graph1()

        ibd_graph = IBDGraphSequence()
        graph_idx = -1

        for ibd_seg in ibd_segs:
            updates_ibd_graph = ibd_graph.add_ibd_segment(ibd_seg)
            if updates_ibd_graph:
                continue

        while ibd_graph.has_incoming_ibd_segments() or ibd_graph.has_active_ibd_segments():
            if ibd_graph.has_at_least_one_edge():
                ibd_graph_conn_comps = ibd_graph.get_connected_components()

                l1 = len(ibd_graph_conn_comps)
                l2 = len(true_conn_comps[graph_idx])
                self.assertTrue(
                    l1 == l2,
                    msg=f"ibd graph {graph_idx} # conn comps {l1} != {l2} # true conn comps"
                )

                for true_comp in true_conn_comps[graph_idx]:
                    true_comp_idx = set(
                        [ibd_graph.get_sample_hap_idx(c[0], c[1]) for c in true_comp]
                    )
                    self.assertTrue(
                        true_comp_idx in ibd_graph_conn_comps,
                        msg=f"missing {true_comp_idx} from graph {graph_idx}"
                    )

            ibd_graph.advance_to_next_breakpoint()
            graph_idx += 1

        self.assertTrue(
            graph_idx == 5,
            msg="failed to iterate through all IBD graphs"
        )


def load_ibd_seg_file(ibd_seg_file):
    ibd_segs = []

    with open(ibd_seg_file) as f:
        for line in f:
            if line[0] == "#": continue
            sid1, hid1, sid2, hid2, chrom, start, end = line.strip("\n").split()
            start = int(start)
            end = int(end)

            seg = IBDSegment(sid1, sid2, hid1, hid2, chrom, start, end)
            ibd_segs.append(seg)
    return ibd_segs


if __name__ == "__main__":
    unittest.main()
