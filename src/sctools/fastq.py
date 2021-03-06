from collections import namedtuple
from . import reader


class Record:
    """Fastq record object

    Defines several properties for accessing fastq record information:
    :property bytes name: name field
    :property bytes sequence: sequence field
    :property bytes name2: second name field
    :property bytes quality: quality field

    Also defines several methods for accessing SEQC annotation fields:
    :property float average_quality: return the mean quality of FastqRecord
    """

    __slots__ = ['_name', '_sequence', '_name2', '_quality']

    def __init__(self, record):
        """
        :param [str|bytes] record: iterable of four strings or bytes objects defining a
          single fastq record
        """
        self.name, self.sequence, self.name2, self.quality = record

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if not isinstance(value, (bytes, str)):
            raise TypeError('fastq name must be str or bytes')
        elif not value.startswith(b'@'):
            raise ValueError('fastq name must start with @')
        else:
            self._name = value

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, value):
        if not isinstance(value, (bytes, str)):
            raise TypeError('fastq sequence must be str or bytes')
        else:
            self._sequence = value

    @property
    def name2(self):
        return self._name2

    @name2.setter
    def name2(self, value):
        if not isinstance(value, (bytes, str)):
            raise TypeError('fastq name2 must be str or bytes')
        else:
            self._name2 = value

    @property
    def quality(self):
        return self._quality

    @quality.setter
    def quality(self, value):
        if not isinstance(value, (bytes, str)):
            raise TypeError('fastq quality must be str or bytes')
        else:
            self._quality = value

    def __bytes__(self):
        return b''.join((self.name, self.sequence, self.name2, self.quality))

    def __str__(self):
        return b''.join((self.name, self.sequence, self.name2, self.quality)).decode()

    def __repr__(self):
        return "Name: %s\nSequence: %s\nName2: %s\nQuality: %s\n" % (
            self.name, self.sequence, self.name2, self.quality)

    def __len__(self):
        return len(self.sequence)

    def average_quality(self):
        """return the average quality of this record"""
        # -33 due to solexa/illumina phred conversion
        return sum(c for c in self.quality[:-1]) / (len(self.quality) - 1) - 33


class StrRecord(Record):
    """Fastq record object

    Defines several properties for accessing fastq record information:
    :property str name: name field
    :property str sequence: sequence field
    :property str name2: second name field
    :property str quality: quality field

    Utilities for parsing read attributes:
    :property average_quality: return the mean quality of FastqRecord
    """

    def __bytes__(self):
        return ''.join((self.name, self.sequence, self.name2, self.quality)).encode()

    def __str__(self):
        return ''.join((self.name, self.sequence, self.name2, self.quality))

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if not isinstance(value, (bytes, str)):
            raise TypeError('fastq name must be str or bytes')
        elif not value.startswith('@'):
            raise ValueError('fastq name must start with @')
        else:
            self._name = value

    def average_quality(self):
        """calculate the average quality of the fastq read

        :return float: average quality of fastq read
        """
        b = self.quality[:-1].encode()
        return sum(c for c in b) / len(b) - 33  # -33 due to solexa/illumina phred conversion



class Reader(reader.Reader):
    """
    Fastq Reader, defines some special methods for reading and summarizing fastq data:

    :method __iter__: Iterator over fastq Record objects
    :method __len__: return number of records in file
    :method estimate_sequence_length: estimate the length of fastq sequences in file
    """

    @staticmethod
    def record_grouper(iterable):
        """
        creates 4 iterators on the same iterable; each moves the pointer forward when called,
        yielding 4 objects at a time
        """
        args = [iter(iterable)] * 4
        return zip(*args)

    def __iter__(self):
        record_type = StrRecord if self._mode == 'r' else Record
        for record in self.record_grouper(super().__iter__()):
            yield record_type(record)


# namedtuple that defines the start and end position of a barcode sequence and provides the name
# for both a quality and sequence tag
EmbeddedBarcode = namedtuple('Tag', ['start', 'end', 'sequence_tag', 'quality_tag'])


class EmbeddedBarcodeGenerator(Reader):

    def __init__(self, embedded_barcodes, *args, **kwargs):
        """
        parse fastq files for barcodes(s) defined by EmbeddedBarcode objects, producing a generator
        that yields extracted barcode objects in a form that is consumable by pysam's bam and sam
        set_tag methods.

        :param [EmbeddedBarcode] embedded_barcodes: list of tag objects defining start and end of
          the sequence containing the tag, plus the string quality and sequence tags
        :param list|str files: (default sys.stdin) fastq file or list of files to be read.
        :param mode: (Default 'r') will return string objects. Change to 'rb' to
          return bytes objects.
        """
        super().__init__(*args, **kwargs)
        self.embedded_barcodes = embedded_barcodes

    def __iter__(self):
        """iterates over barcodes extracted from fastq"""
        for record in super().__iter__():  # iterates records; we extract barcodes.
            barcodes = []
            for barcode in self.embedded_barcodes:
                barcodes.extend(self.extract_barcode(record, barcode))
            yield barcodes

    @staticmethod
    def extract_barcode(record, eb):
        """extracts barcodes from a fastq record according to an EmbeddedBarcode object.

        :param FastqRecord record: record to extract from
        :param EmbeddedBarcode eb: defines the barcode start and end and the names of the quality
          and sequence tags
        :return tuple: (sequence tag name, tag content (nucleotides), tag type (Z))
        :return tuple: (quality tag name, tag content (qualities), tag type (Z))
        """
        seq = record.sequence[eb.start:eb.end]
        qual = record.quality[eb.start:eb.end]
        return (eb.sequence_tag, seq, 'Z'), (eb.quality_tag, qual, 'Z')
