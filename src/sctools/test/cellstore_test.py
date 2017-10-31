import unittest
from io import BytesIO
import numpy as np
from google.cloud import storage
from .. import cellstore

# todo fix these tests


class TestCellStore(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.rs = cellstore.CellStore._from_test(1000)
        cls.archive_prefix = 'read_store_test'

    def test_setup(self):
        self.assertTrue(True)  # just testing class setup

    def test_diff_recordarray(self):
        self.rs.index('cell')

    def test_create_save_index(self):
        self.rs.save(self.archive_prefix)

    def test_get_header(self):
        h = self.rs.get_header(self.archive_prefix + '.npi')
        print(h)

    def test_figure_out_numpy_save(self):
        from numpy.lib import format
        print('serialized dtype')
        print(format.dtype_to_descr(self.rs.data.dtype))

        print('header data')
        d = format.header_data_from_array_1_0(self.rs.data)
        print(d)

        print('header magic string')
        fp = BytesIO()
        format._write_array_header(fp, d)
        fp.seek(0); print(fp.read())

        with open(self.archive_prefix + '.npy', 'rb') as f:
            print(format.read_array_header_1_0(f))

    def test_load(self):
        self.rs.sort('cell')
        rs = cellstore.CellStore.load(self.archive_prefix + '.npy')
        print(self.rs.data.dtype, self.rs.data.shape)
        print(rs.data.dtype, rs.data.shape)
        print(rs.data[:10])
        print(self.rs.data[:10])
        self.assertTrue(np.array_equal(rs.data, self.rs.data))

    def test_load_range(self):
        self.rs.save(self.archive_prefix)
        range_ = slice(20, 25)
        rs2 = cellstore.CellStore.from_range(self.archive_prefix + '.npy', range_)
        print(rs2.data)
        print(self.rs.data[20:25])
        self.assertTrue(np.array_equal(self.rs.data[range_], rs2))

    def test_load_range_gcloud(self):
        filename = 'gs://broad-dsde-mint-dev-teststorage/readstore/read_store_test.npy'
        range_ = slice(20, 25)
        client = storage.Client(project='broad-dsde-mint-dev')  # requires previous OAuth login
        print(cellstore.CellStore.from_gstorage_bucket(filename, range_, client))

    def test_blosc_compression_ratio(self):
        rs = cellstore.CellStore._from_test(int(1e7))  # 1m * 6 bytes
        print(rs.data.nbytes)
        result = rs.compress_blosc()
        print(len(result))
        # ~ 3x in worst case scenario (random, noise). Not bad!

    @classmethod
    def tearDownClass(cls):
        pass  # remove any saved files


if __name__ == "__main__":
    unittest.main()
