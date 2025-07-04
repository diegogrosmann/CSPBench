from csp_blfga.utils.distance import hamming_distance, max_hamming


def test_hamming_distance():
    assert hamming_distance("AAAA", "AAAT") == 1
    assert hamming_distance("AAAA", "TTTT") == 4
    assert hamming_distance("ATAT", "TATA") == 4


def test_max_hamming():
    strings = ["AAAA", "AAAT", "AATT", "TTTT"]
    assert max_hamming("AAAA", strings) == 4
    assert max_hamming("TTTT", strings) == 4
