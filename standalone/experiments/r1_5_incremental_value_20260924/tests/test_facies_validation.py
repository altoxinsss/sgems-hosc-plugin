from pathlib import Path

import numpy as np

from validation import channel_count_2x2x2, channel_facies_mask, downsample_facies_mode, read_facies_dat, strict_majority_channel_labels, subset_facies_by_grid_coordinates


def test_facies_reader_and_blockwise_mode_downsampling(tmp_path: Path):
    # Each 2x2x2 block has a categorical majority; the result must never be an
    # arithmetic mean of facies identifiers.
    fine = np.array(
        [
            [
                [2, 2],
                [2, 1],
                [1, 1],
                [1, 2],
            ],
            [
                [2, 2],
                [3, 2],
                [1, 3],
                [1, 1],
            ],
        ],
        dtype=np.int16,
    )
    path = tmp_path / "facies.dat"
    path.write_text("title\n1\nfacies\n" + " ".join(str(value) for value in fine.reshape(-1)) + "\n", encoding="utf-8")
    loaded = read_facies_dat(path, fine.shape)
    assert np.array_equal(loaded, fine)
    coarse = downsample_facies_mode(loaded, (1, 2, 1))
    assert np.array_equal(coarse, np.array([[[2], [1]]], dtype=np.int16))
    assert np.array_equal(channel_facies_mask(coarse), np.array([[[True], [False]]]))


def test_one_based_z62_coordinate_selects_numpy_index_61():
    facies = np.arange(100, dtype=np.int16).reshape(100, 1, 1)
    subset = subset_facies_by_grid_coordinates(facies, np.array([1]), np.array([1]), np.array([62]))
    assert subset.shape == (1, 1, 1)
    assert subset[0, 0, 0] == 61


def test_strict_channel_majority_excludes_exact_ties():
    facies = np.zeros((2, 2, 4), dtype=np.int16)
    facies[..., :2] = 2
    facies[0, 0, 2:] = 2
    facies[0, 1, 2:] = 2
    counts = channel_count_2x2x2(facies)
    labels, valid = strict_majority_channel_labels(counts)
    assert np.array_equal(counts, np.array([[[8, 4]]], dtype=np.int8))
    assert np.array_equal(labels, np.array([[[True, False]]]))
    assert np.array_equal(valid, np.array([[[True, False]]]))
