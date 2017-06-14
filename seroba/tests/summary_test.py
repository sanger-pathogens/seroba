import unittest
import os
import filecmp
import shutil
import tempfile
from seroba import serotyping

modules_dir = os.path.dirname(os.path.abspath(serotyping.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestSummary(unittest.TestCase):
    def test__summary(self):
        pass
