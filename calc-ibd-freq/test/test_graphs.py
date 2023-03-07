import numpy as np

# IBD seg file
# #SID1   HID2    SID2    HID2    CHROM   START   END
# 1       1       2       1       1       1       10
# 1       1       4       1       1       1       10
# 2       1       3       1       1       5       15
# 1       1       3       2       1       20      25

def get_test_graph1():
    """Test graphs corresponding to test/data/ibd_segs.txt
    """


    true_conn_comps = [
        # from position 1-4
        [
            [("1", "1"), ("2", "1"), ("4", "1")]
        ],
        # from position 5-10
        [
            [("1", "1"), ("2", "1"), ("4", "1"), ("3", "1")]
        ],
        # from position 11-15
        [
            [("2", "1"), ("3", "1")]
        ],
        # from 16-19
        [
            []
        ],
        # from position 20-25
        [
            [("1", "1"), ("3", "2")]
        ],
        # from 26-end
        [
            []
        ]
    ]
    true_spans = [
        [1, 4], [5, 10], [11, 15], [16, 19], [20, 25], [26, np.inf]
    ]

    return true_conn_comps, true_spans

# IBD seg file
# #SID1   HID2    SID2    HID2    CHROM   START   END
# 1       1       2       1       1       5       10
# 1       2       2       2       1       5       10
def get_test_graph2():
    """Test graphs corresponding to test/data/ibd_segs_2.txt
    """

    true_conn_comps = [
        [],
        [
            [("1", "1"), ("2", "1")],
            [("1", "2"), ("2", "2")]
        ],
        []
    ]

    true_spans = [
        [1, 4],
        [5, 10],
        [11, np.inf]
    ]

    return true_conn_comps, true_spans

