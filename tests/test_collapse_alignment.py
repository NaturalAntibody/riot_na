from riot_na.data.model import AlignmentString
from riot_na.schemes.collapse_alignment import (
    _collapse_ins_del_ordered,
    collapse_alignment_str,
    collapse_ins_del,
)


def test_collapse_alignment_str():
    # given
    alignment_str = AlignmentString("MDDIM")
    # when
    result = collapse_alignment_str(alignment_str)

    # then
    assert result == "MMDM"

    # given
    alignment_str = AlignmentString("MDDDIIIIIIM")
    # when
    result = collapse_alignment_str(alignment_str)

    # then
    assert result == "MMMIIIMM"

    # given
    alignment_str = AlignmentString("MDDDMIDIDIIM")
    # when
    result = collapse_alignment_str(alignment_str)

    # then
    assert result == "MDDDMMIIMM"

    # indels in the middle
    assert collapse_alignment_str(AlignmentString("MMMIIIIDDDMMM")) == "MMMMMIMMMM"

    # indels ordered
    assert collapse_alignment_str(AlignmentString("MMMIIIIDDDMMM"), ordered=True) == "MMMIMMMMMM"


def test_collapse_ins_del_ordered():
    assert _collapse_ins_del_ordered(list("DDDIIII")) == "MMMI"
    assert _collapse_ins_del_ordered(list("IIIDDDD")) == "MMMD"

    assert _collapse_ins_del_ordered(list("DDDDIII")) == "DMMM"
    assert _collapse_ins_del_ordered(list("IIIIDDD")) == "IMMM"

    assert _collapse_ins_del_ordered(list("IDDDI")) == "MDM"
    assert _collapse_ins_del_ordered(list("DIIID")) == "MIM"

    assert _collapse_ins_del_ordered(list("DDDIII")) == "MMM"
    assert _collapse_ins_del_ordered(list("IIIDDD")) == "MMM"

    assert _collapse_ins_del_ordered(list("IDIDID")) == "MMM"

    assert _collapse_ins_del_ordered(list("IDIIDIDDI")) == "MIMMM"

    assert _collapse_ins_del_ordered(list("IIIIDDDIIDIDDDIID")) == "IMMMIMMMMD"

    assert _collapse_ins_del_ordered(list("DDIDIDID")) == "DMMMD"

    assert _collapse_ins_del_ordered(list("IIDIDIDI")) == "IMMMI"


def test_collapse_insdel_insertions_in_junctions():
    assert collapse_ins_del(list("IDDDDDDDD ")) == "MDDDDDDD"

    assert collapse_ins_del(list("DIDDDDDDD ")) == "MDDDDDDD"

    assert collapse_ins_del(list("DIDDDDDID ")) == "MDDDDDM"

    assert collapse_ins_del(list("IDDDDDDID ")) == "MDDDDDM"

    assert collapse_ins_del(list("IDDDDDDDI ")) == "MDDDDDM"

    assert _collapse_ins_del_ordered(list("DIDDDDDID")) == "MDDDDMD"

    assert _collapse_ins_del_ordered(list("IMDDDDDDD")) == "MMDDDDDD"
    assert _collapse_ins_del_ordered(list("IMDDDDDMI")) == "MMDDDMM"


if __name__ == "__main__":
    test_collapse_alignment_str()
    test_collapse_ins_del_ordered()
    test_collapse_insdel_insertions_in_junctions()
