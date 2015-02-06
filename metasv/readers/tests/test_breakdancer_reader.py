from ..breakdancer_reader import BreakDancerRecord
import unittest

class TestBreakDancerRecord(unittest.TestCase):

    def setUp(self):
        self.breakdancer_string = "22   26891493	4+0-	22	26891659	0+4-	DEL	224	88	4	jobs/617dfc0a-f87f-4a74-83d3-ebb72bf3659d/15.recalibrated.bam|4	NA	NA	NA	NA	NA	NA	NA	1.41	NA	NA	NA	NA	NA	NA	NA	NA"

    def test_parsing(self):
        record = BreakDancerRecord(self.breakdancer_string)

        self.assertEqual(record.sv_len, 224)


if __name__ == '__main__':
    unittest.main()

