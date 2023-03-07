import numpy as np
import unittest

from src.ibd_freq_estimator.ibd_file import IBDFile
from test.test_graphs import get_test_graph1, get_test_graph2

class TestIBDFile(unittest.TestCase):

    def test_ibd_file(self):
        ibd_files = [
            "./test/data/ibd_segs.txt",
            "./test/data/ibd_segs_2.txt"
        ]
        test_graphs = [
            get_test_graph1(),
            get_test_graph2()
        ]

        for ibd_file_path, test_graph in zip(ibd_files, test_graphs):
            ibd_file = IBDFile(
                ibd_file_path,
                sid1_col=0,
                sid2_col=2,
                hid1_col=1,
                hid2_col=3,
                chr_col=4,
                start_pos_col=5,
                end_pos_col=6,
                has_header=True
            )

            true_conn_comps, true_spans = test_graph
            for i,ibd_graph in enumerate(ibd_file):
                span = ibd_graph.get_span()
                conn_comps = ibd_graph.get_connected_components()

                for true_comp in true_conn_comps[i]:
                    true_comp_idx = set(
                        [ibd_graph.get_sample_hap_idx(c[0], c[1]) for c in true_comp]
                    )

                    self.assertTrue(
                        true_comp_idx == set() or true_comp_idx in conn_comps,
                        msg=f"missing {true_comp_idx} from graph {i} ({ibd_file_path})"
                    )

                    self.assertTrue(
                        true_spans[i][0] == span[0],
                        msg=f"incorrect start span {true_spans[i][0]} != {span[0]} ({ibd_file_path})"
                    )

                    self.assertTrue(
                        true_spans[i][1] == span[1],
                        msg=f"incorrect end span {true_spans[i][1]} != {span[1]} ({ibd_file_path})"
                    )


if __name__ == "__main__":
    unittest.main()
