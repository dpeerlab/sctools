import numpy as np
import pickle
from io import BytesIO
from google.resumable_media.requests import Download
import blosc

# define a trivial data structure that:
# 1. slurps from .bam
# 2. is indexed
# 3. through the index, is randomly accessible


# minimal specification for read 10x data
_10x_dtype = [
    ('cell', np.uint32),  # 4 bytes
    ('cell_quality', np.dtype('(16,)b')),  # 16 bytes
    ('rmt', np.uint16),  # 2 bytes
    ('rmt_quality', np.dtype('(8,)b')),  # 8 bytes
    ('sequence1', np.dtype('(25,)b')),  # 25 bytes (24 + change)
    ('quality', np.dtype('(98,)b')),  # 98 bytes
    ('flag', np.uint8),  # 1 byte
    ('chromosome', np.uint8),  # 1 byte
    ('strand', np.uint8),  # 1 byte
    ('position', np.uint32),  # 4 bytes
    # ('CIGAR', np.dtype('(???)')),  # ??? how many to give this?
]

# assuming 200m reads, sum = 2.6GB for:
# cell, rmt, flag, chrom, strand, position

# sum = ~ 32 gb for the whole thing above, but compressible down to


class CellStore:

    # test specification
    _dtype = [
        ('cell', np.uint32),
        ('umi', np.uint16)
    ]
    _record_size = (32 + 16) // 8  # todo automate this!

    def __init__(self, data):
        self._data = data

    @property
    def data(self):
        return self._data

    @classmethod
    def from_bam(cls, bamfile):
        """

        :param bamfile:
        :return:
        """
        raise NotImplementedError

    @classmethod
    def _from_test(cls, n_records):
        data = np.zeros(n_records, cls._dtype)
        data['cell'] = np.random.randint(0, 32, n_records, dtype=cls._dtype[0][1])
        data['umi'] = np.random.randint(0, 16, n_records, dtype=cls._dtype[1][1])
        return cls(data)

    @classmethod
    def from_range(cls, archive, range_):
        """

        :param archive:
        :param range_:
        :return:
        """
        # todo this is hardcoded right now for testing simplicity
        # need to know header length of original file (this should be static, and based on this
        # archive.
        # magic len = 6 (prefix) + 2
        # header length identifier (2)
        # header length (102)
        # for this data, an offset of 112 should get us into the file.
        with open(archive, 'rb') as fp:
            fp.seek(112 + range_.start * cls._record_size)
            array = np.fromfile(fp, dtype=cls._dtype, count=range_.stop - range_.start)
        return cls(array)

    def sort(self, field):
        self.data.sort(kind='mergesort', order=field)

    def index(self, field):
        """simple range index over grouped objects using a pandas series

        :param field:
        :return:
        """

        self.data.sort(kind='mergesort', order=field)
        flag = np.concatenate(([True], self.data[field][1:] != self.data[field][:-1]))
        values = self.data[field][flag]

        # identify the locations where changes between values occur, add end of array
        index = np.concatenate([np.where(flag)[0], [self.data.shape[0]]])

        # create slice objects for each range
        ranges = [slice(i, j) for i, j in zip(index[:-1][:, None], index[1:][:, None])]

        # could be changed to np.array and stored in .npz file
        return dict(zip(values, ranges))

    def save(self, filestem):
        cell_index = self.index('cell')
        with open(filestem + '.npi', 'wb') as f:
            pickle.dump(cell_index, f)
        np.save(filestem + '.npy', self.data)

    @staticmethod
    def get_header(archive):
        with open(archive, 'rb') as f:
            data = f.readlines()
        print(data[0:20])

    @classmethod
    def load(cls, archive):
        data = np.load(archive)
        return cls(data)

    @classmethod
    def from_gstorage_bucket(cls, filename, range_, authenticated_client):

        # set up bucket and blob
        prefix, _, bucket, *blobname = filename.split('/')
        blobname = '/'.join(blobname)
        bucket = authenticated_client.bucket(bucket)
        blob = bucket.blob(blobname)

        # create inputs, determine byte range
        file_obj = BytesIO()
        start_byte = 112 + range_.start * cls._record_size
        stop_byte = 112 + range_.stop * cls._record_size - 1

        # use google's internals to format the request, but stick our header in
        download_url = blob._get_download_url()
        transport = blob._get_transport(authenticated_client)
        download = Download(
            download_url,
            stream=file_obj,
            headers={'range': 'bytes={s}-{e}'.format(s=start_byte, e=stop_byte)})
        download.consume(transport)

        # could chunk this with ChunkedDownload instead of Download
        # and code from numpy.load (they use chunking to read)
        file_obj.seek(0)
        array = np.frombuffer(file_obj.read(), dtype=cls._dtype)

        return cls(array)

    def subset_from_index(self):
        raise NotImplementedError  # should be trivial to do this!

    def compress_blosc(self):
        return blosc.compress_ptr(
            self.data.__array_interface__['data'][0],
            self.data.size,
            self.data.dtype.itemsize, 5, True)
