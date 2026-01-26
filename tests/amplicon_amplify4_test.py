# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Amplify 4-specific tests for AmpliconGenerator."""

import pytest

from amplifyp.amplicon import AmpliconGenerator
from amplifyp.dna import DNADirection
from amplifyp.repliconf import Repliconf
from amplifyp.settings import DEFAULT_REPLICATION_SETTINGS
from tests.examples.amplify4_examples import (
    amplify4_circular_example,
    amplify4_linear_example,
    fragment_177bp,
    fragment_190bp,
    fragment_220bp,
    fragment_566bp,
    fragment_660bp,
    fragment_1250bp,
    fragment_1280bp,
    fragment_1626bp,
    fragment_2033bp,
    fragment_2516bp,
    primer_11bp,
    primer_1701,
    primer_10289,
    primer_10290,
)


def test_amplify4_11bp_1701_10289_10290() -> None:
    """Test amplicon generation using Amplify 4 examples."""
    amplicon_generator = AmpliconGenerator(amplify4_linear_example)

    repliconf_11bp = Repliconf(
        amplify4_linear_example, primer_11bp, DEFAULT_REPLICATION_SETTINGS
    )
    repliconf_1701 = Repliconf(
        amplify4_linear_example, primer_1701, DEFAULT_REPLICATION_SETTINGS
    )
    repliconf_10289 = Repliconf(
        amplify4_linear_example, primer_10289, DEFAULT_REPLICATION_SETTINGS
    )
    repliconf_10290 = Repliconf(
        amplify4_linear_example, primer_10290, DEFAULT_REPLICATION_SETTINGS
    )

    amplicon_generator.add(repliconf_11bp)
    amplicon_generator.add(repliconf_1701)
    amplicon_generator.add(repliconf_10289)
    amplicon_generator.add(repliconf_10290)

    amplicons = amplicon_generator.get_amplicons()
    assert len(amplicons) == 7

    encountered_177 = False
    encountered_190 = False
    encountered_220 = False
    encountered_566 = False
    encountered_660 = False
    encountered_2033 = False
    encountered_2516 = False

    for amplicon in amplicons:
        if len(amplicon.product.seq) == 177:
            assert amplicon.product.seq.upper() == fragment_177bp.seq.upper()
            assert amplicon.q_score == 146
            encountered_177 = True
        elif len(amplicon.product.seq) == 190:
            assert amplicon.product.seq.upper() == fragment_190bp.seq.upper()
            # N.B. This Q-score is different than what is reported by Amplify 4.
            assert amplicon.q_score == pytest.approx(699.9337337142201)
            encountered_190 = True
        elif len(amplicon.product.seq) == 220:
            assert amplicon.product.seq.upper() == fragment_220bp.seq.upper()
            assert amplicon.q_score == 189
            encountered_220 = True
        elif len(amplicon.product.seq) == 566:
            assert amplicon.product.seq.upper() == fragment_566bp.seq.upper()
            # N.B. This Q-score is different than what is reported by Amplify 4.
            # See test_10290_10289_566bp_manual for more information.
            assert amplicon.q_score == pytest.approx(55696.52822187758)
            encountered_566 = True
        elif len(amplicon.product.seq) == 660:
            assert amplicon.product.seq.upper() == fragment_660bp.seq.upper()
            assert amplicon.q_score == 629
            encountered_660 = True
        elif len(amplicon.product.seq) == 2033:
            assert amplicon.product.seq.upper() == fragment_2033bp.seq.upper()
            assert amplicon.q_score == 1993
            encountered_2033 = True
        elif len(amplicon.product.seq) == 2516:
            assert amplicon.product.seq.upper() == fragment_2516bp.seq.upper()
            assert amplicon.q_score == 2476
            encountered_2516 = True
        else:
            raise AssertionError(
                f"Unexpected amplicon length: {len(amplicon.product.seq)}"
            )

    assert encountered_177
    assert encountered_190
    assert encountered_220
    assert encountered_566
    assert encountered_660
    assert encountered_2033
    assert encountered_2516


def test_10290_10289_566bp_manual() -> None:
    """Verify 10290_10289_566bp manually.

    This calculates the 566bp amplicon generated from the built-in template
    using the built-in primers with the label 10290 and 10289. The target
    sequence length is 526bp. The whole amplicon is 566bp.

    The Q-score here is different than what is reported by Amplify 4. This
    is because in Amplify 4, the replication origin quality score is calculated
    from truncated primability and stability scores. In this implementation,
    the replication origin quality score is calculated from the full-length
    primability and stability scores.
    """
    repliconf_10289 = Repliconf(
        amplify4_linear_example, primer_10289, DEFAULT_REPLICATION_SETTINGS
    )
    repliconf_10290 = Repliconf(
        amplify4_linear_example, primer_10290, DEFAULT_REPLICATION_SETTINGS
    )
    repliconf_10289.search()
    repliconf_10290.search()
    assert repliconf_10290.origin_db.fwd == [710, 2177]
    origin_fwd = repliconf_10290.origin_from_db(DNADirection.FWD, 1)
    assert origin_fwd.primability == pytest.approx(0.8059701492537313)
    assert origin_fwd.stability == pytest.approx(0.4717741935483871)
    assert origin_fwd.quality == pytest.approx(0.09718042850264788)
    origin_rev = repliconf_10289.origin_from_db(DNADirection.REV, 0)
    assert origin_rev.primability == pytest.approx(1.0)
    assert origin_rev.stability == pytest.approx(1.0)
    assert origin_rev.quality == pytest.approx(1.0)
    fwd_quality = origin_fwd.quality
    rev_quality = origin_rev.quality
    assert 526 / (fwd_quality * rev_quality) ** 2 == pytest.approx(
        55696.52822187758
    )


def test_amplify4_circular_11bp_1701_10289_10290() -> None:
    """Test amplicon generation using Amplify 4 examples."""
    amplicon_generator = AmpliconGenerator(amplify4_circular_example)

    repliconf_11bp = Repliconf(
        amplify4_circular_example, primer_11bp, DEFAULT_REPLICATION_SETTINGS
    )
    repliconf_1701 = Repliconf(
        amplify4_circular_example, primer_1701, DEFAULT_REPLICATION_SETTINGS
    )
    repliconf_10289 = Repliconf(
        amplify4_circular_example, primer_10289, DEFAULT_REPLICATION_SETTINGS
    )
    repliconf_10290 = Repliconf(
        amplify4_circular_example, primer_10290, DEFAULT_REPLICATION_SETTINGS
    )

    amplicon_generator.add(repliconf_11bp)
    amplicon_generator.add(repliconf_1701)
    amplicon_generator.add(repliconf_10289)
    amplicon_generator.add(repliconf_10290)

    amplicons = amplicon_generator.get_amplicons()
    assert len(amplicons) == 10

    encountered_177 = False
    encountered_190 = False
    encountered_220 = False
    encountered_566 = False
    encountered_660 = False
    encountered_2033 = False
    encountered_2516 = False
    # Circular amplicons
    encountered_1250 = False
    encountered_1280 = False
    encountered_1626 = False

    for amplicon in amplicons:
        if len(amplicon.product.seq) == 177:
            assert amplicon.product.seq.upper() == fragment_177bp.seq.upper()
            assert amplicon.q_score == 146
            encountered_177 = True
        elif len(amplicon.product.seq) == 190:
            assert amplicon.product.seq.upper() == fragment_190bp.seq.upper()
            # N.B. This Q-score is different than what is reported by Amplify 4.
            assert amplicon.q_score == pytest.approx(699.9337337142201)
            encountered_190 = True
        elif len(amplicon.product.seq) == 220:
            assert amplicon.product.seq.upper() == fragment_220bp.seq.upper()
            assert amplicon.q_score == 189
            encountered_220 = True
        elif len(amplicon.product.seq) == 566:
            assert amplicon.product.seq.upper() == fragment_566bp.seq.upper()
            # N.B. This Q-score is different than what is reported by Amplify 4.
            # See test_10290_10289_566bp_manual for more information.
            assert amplicon.q_score == pytest.approx(55696.52822187758)
            encountered_566 = True
        elif len(amplicon.product.seq) == 660:
            assert amplicon.product.seq.upper() == fragment_660bp.seq.upper()
            assert amplicon.q_score == 629
            encountered_660 = True
        elif len(amplicon.product.seq) == 2033:
            assert amplicon.product.seq.upper() == fragment_2033bp.seq.upper()
            assert amplicon.q_score == 1993
            encountered_2033 = True
        elif len(amplicon.product.seq) == 2516:
            assert amplicon.product.seq.upper() == fragment_2516bp.seq.upper()
            assert amplicon.q_score == 2476
            encountered_2516 = True
        # Circular amplicons
        elif len(amplicon.product.seq) == 1250:
            assert amplicon.product.seq.upper() == fragment_1250bp.seq.upper()
            # N.B. This Q-score is different than what is reported by Amplify 4.
            assert amplicon.q_score == pytest.approx(5405.777515729951)
            encountered_1250 = True
        elif len(amplicon.product.seq) == 1280:
            assert amplicon.product.seq.upper() == fragment_1280bp.seq.upper()
            assert amplicon.q_score == 1258
            encountered_1280 = True
        elif len(amplicon.product.seq) == 1626:
            assert amplicon.product.seq.upper() == fragment_1626bp.seq.upper()
            # N.B. This Q-score is different than what is reported by Amplify 4.
            assert amplicon.q_score == pytest.approx(168889.6625739444)
            encountered_1626 = True
        else:
            raise AssertionError(
                f"Unexpected amplicon length: {len(amplicon.product.seq)}"
            )

    assert encountered_177
    assert encountered_190
    assert encountered_220
    assert encountered_566
    assert encountered_660
    assert encountered_2033
    assert encountered_2516
    # Circular amplicons
    assert encountered_1250
    assert encountered_1280
    assert encountered_1626
