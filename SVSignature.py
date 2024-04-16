import logging

class Signature:
    """Signature class for basic signatures of structural variants. An signature is always detected from a single read.
    """

    def __init__(self, contig, start, end, signature, read):
        self.contig = contig
        self.start = start
        self.end = end
        self.signature = signature
        self.read = read
        self.type = None
        self.data = None
        self.svlen = self.end - self.start
        if self.end < self.start:
            logging.warning("Signature with invalid coordinates (end < start): " + self.as_string())

    def get_source(self):
        return (self.contig, self.start, self.end)

    def get_key(self):
        contig, start, end = self.get_source()
        return (self.type, contig, end)

    def downstream_distance_to(self, signature2):
        """Return distance >= 0 between this signature's end and the start of signature2."""
        this_contig, this_start, this_end = self.get_source()
        other_contig, other_start, other_end = signature2.get_source()
        if self.type == signature2.type and this_contig == other_contig:
            return max(0, other_start - this_end)
        else:
            return float("inf")

class SignatureDeletion(Signature):
    """SV Signature: a region (contig:start-end) has been deleted and is not present in sample"""

    def __init__(self, contig, start, end, signature, read):
        self.contig = contig
        assert end >= start
        # 0-based start of the deletion (first deleted base)
        self.start = start
        # 0-based end of the deletion (one past the last deleted base)
        self.end = end
        self.signature = signature
        self.read = read
        self.type = "DEL"
        self.svlen = self.end - self.start

class SignatureInsertion(Signature):
    """SV Signature: a region of length end-start has been inserted at contig:start"""

    def __init__(self, contig, start, end, signature, read):
        self.contig = contig
        assert end >= start
        # 0-based start of the insertion (base after the insertion)
        self.start = start
        # 0-based start of the insertion (base after the insertion) + length of the insertion
        self.end = end
        self.signature = signature
        self.read = read
        self.type = "INS"
        self.svlen = self.end - self.start

    def get_key(self):
        contig, start, end = self.get_source()
        return (self.type, contig, start)

    def downstream_distance_to(self, signature2):
        """Return distance >= 0 between this signature's end and the start of signature2."""
        this_contig, this_start, this_end = self.get_source()
        other_contig, other_start, other_end = signature2.get_source()
        if self.type == signature2.type and this_contig == other_contig:
            return max(0, other_start - this_start)
        else:
            return float("inf")

class SignatureInversion(Signature):
    """SV Signature: a region (contig:start-end) has been inverted in the sample"""

    def __init__(self, contig, start, end, signature, read, direction):
        self.contig = contig
        assert end >= start
        # 0-based start of the inversion (first inverted base)
        self.start = start
        # 0-based end of the inversion (one past the last inverted base)
        self.end = end
        self.signature = signature
        self.read = read
        self.type = "INV"
        self.direction = direction
        self.svlen = self.end - self.start

class SignatureDuplicationTandem(Signature):
    """SV Signature: a region (contig:start-end) has been tandemly duplicated"""

    def __init__(self, contig, start, end, signature, read):
        self.contig = contig
        assert end >= start
        # 0-based start of the region (first copied base)
        self.start = start
        # 0-based end of the region (one past the last copied base)
        self.end = end
        self.signature = signature
        self.read = read
        self.type = "DUP"
        self.svlen = self.end - self.start

class SignatureTranslocation(Signature):
    """SV Signature: two positions (contig1:pos1 and contig2:pos2) are connected in the sample"""

    def __init__(self, contig1, pos1, direction1, contig2, pos2, direction2, signature, read):
        if contig1 < contig2 or (contig1 == contig2 and pos1 < pos2):
            self.contig1 = contig1
            # 0-based source of the translocation (first base before the translocation)
            self.pos1 = pos1
            self.source_direction = direction1
            self.contig2 = contig2
            # 0-based destination of the translocation (first base after the translocation)
            self.pos2 = pos2
            self.dest_direction = direction2
        else:
            print('translocation error')
        self.signature = signature
        self.read = read
        self.type = "BND"
        self.svlen = None

    def get_source(self):
        return (self.contig1, self.pos1, self.pos1 + 1)

    def get_destination(self):
        return (self.contig2, self.pos2, self.pos2 + 1)

    def get_key(self):
        return (self.type, self.contig1, self.pos1)



