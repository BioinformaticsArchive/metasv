from abc import ABCMeta, abstractmethod

class NativeSVRecord():
    __metaclass__ = ABCMeta

    @abstractmethod
    def to_vcf_record(self, sample=None):
        pass

    @abstractmethod
    def to_sv_interval(self):
        pass