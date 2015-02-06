import logging
import vcf
from native_sv_reader import NativeSVRecord
from sv_interval import SVInterval

logger = logging.getLogger(__name__)
logging.basicConfig()

class BreakDancerRecord(NativeSVRecord):
    name = "BreakDancer"
    source = set(name)
    valid_svs = set(["DEL", "INS", "INV"])

    def __init__(self, record_string):
        fields = record_string.split()
        self.chr1 = fields[0]
        self.pos1 = int(fields[1])
        self.ori1 = fields[2]
        self.chr2 = fields[3]
        self.pos2 = int(fields[4])
        self.ori2 = fields[5]
        self.sv_type = fields[6]
        self.sv_len = abs(int(fields[7]))
        self.score = float(fields[8])
        self.supporting_read_pairs = int(fields[9])
        self.supporting_reads_pairs_lib = dict(
            map(lambda l: (l[0], int(l[1])), (s.split("|") for s in fields[10].split(":"))))
        self.info = {
            "BD_CHR1": self.chr1,
            "BD_POS1": self.pos1,
            "BD_ORI1": self.ori1,
            "BD_CHR2": self.chr2,
            "BD_POS2": self.pos2,
            "BD_ORI2": self.ori2,
            "BD_SCORE": self.score,
            "BD_SUPPORTING_READ_PAIRS": self.supporting_read_pairs
        }

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return "<" + self.__class__.__name__ + " " + str(self.__dict__) + ">"

    def to_sv_interval(self):
        if self.sv_type not in self.valid_svs:
            logger.error("Unsupported SV type %s" % self.sv_type)
            return None

        cipos = [0, self.pos2 - self.pos1]
        if self.sv_type != "INS": cipos[1] -= abs(self.sv_len)

        return SVInterval(self.chr1,
                          self.pos1 + 1,
                          self.pos2,
                          name=self.name,
                          sv_type=self.sv_type,
                          length=self.sv_len,
                          sources=self.source,
                          cipos=cipos,
                          info=self.info,
                          native_sv=self)

    def to_vcf_record(self, sample=None):
        alt = ["<%s>" % self.sv_type]
        sv_len = -self.sv_len if self.sv_type == "DEL" else self.sv_len
        info = {"SVLEN": sv_len,
                "SVTYPE": self.sv_type}
        if self.sv_type == "DEL" or self.sv_type == "INV":
            info["END"] = self.pos1 + self.sv_len
        elif self.sv_type == "INS":
            info["END"] = self.pos1
        else:
            return None

        info.update(self.info)

        sample_indexes = [0] if sample else []
        sample_calls = [vcf.model._Call(None, sample, vcf.model.make_calldata_tuple("GT")(GT="1/1"))] if sample else []
        return vcf.model._Record(self.chr1,
                                 self.pos1,
                                 ".",
                                 "N",
                                 alt,
                                 ".",
                                 ".",
                                 info,
                                 "GT",
                                 sample_indexes,
                                 sample_calls)

