from csp_blfga.utils.distance import hamming_distance, max_distance


def test_hamming_distance():
    assert hamming_distance("AAAA", "AAAT") == 1
    assert hamming_distance("AAAA", "TTTT") == 4
    assert hamming_distance("ATAT", "TATA") == 4


def test_max_distance():
    strings = ["AAAA", "AAAT", "AATT", "TTTT"]
    assert max_distance("AAAA", strings) == 4
    assert max_distance("TTTT", strings) == 4
