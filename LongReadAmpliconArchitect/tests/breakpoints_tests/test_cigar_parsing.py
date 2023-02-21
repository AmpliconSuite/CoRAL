"""
A module for performing tests on breakpoint_graph/cigar_parsing.py
"""
import unittest


from long_read_aa.breakpoints import cigar_parsing


class TestCigarParsing(unittest.TestCase):

	def test_cigar2posSM(self):
		# Test the function cigar2posSM
		c1, strand1, rl1 = "856S542M", '-', 1398
		c2, strand2, rl2 = "215S202M", '+', 417
        
		qs1, qe1, al1 = cigar_parsing.cigar2posSM(c1, strand1, rl1)
		qs2, qe2, al2 = cigar_parsing.cigar2posSM(c2, strand2, rl2)

		# Expected output
		eqs1, eqe1, eal1 = 0, 541, 542
		eqs2, eqe2, eal2 = 215, 416, 202

		self.assertEqual(qs1, eqs1)
		self.assertEqual(qe1, eqe1)
		self.assertEqual(al1, eal1)
		self.assertEqual(qs2, eqs2)
		self.assertEqual(qe2, eqe2)
		self.assertEqual(al2, eal2)


	def test_cigar2posMS(self):
		# Test the function cigar2posMS
		c3, strand3, rl3 = "437M843S", '-', 1280 
		c4, strand4, rl4 = "980M1770S", '+', 2750
        
		qs3, qe3, al3 = cigar_parsing.cigar2posMS(c3, strand3, rl3)
		qs4, qe4, al4 = cigar_parsing.cigar2posMS(c4, strand4, rl4)

		# Expected output
		eqs3, eqe3, eal3 = 843, 1279, 437
		eqs4, eqe4, eal4 = 0, 979, 980

		self.assertEqual(qs3, eqs3)
		self.assertEqual(qe3, eqe3)
		self.assertEqual(al3, eal3)
		self.assertEqual(qs4, eqs4)
		self.assertEqual(qe4, eqe4)
		self.assertEqual(al4, eal4)


	def test_cigar2posSMS(self):
		# Test the function cigar2posSMS
		c5, strand5, rl5 = "2257S377M42S", '-', 2676 
		c6, strand6, rl6 = "2647S1079M3351S", '+', 7077 
        
		qs5, qe5, al5 = cigar_parsing.cigar2posSMS(c5, strand5, rl5)
		qs6, qe6, al6 = cigar_parsing.cigar2posSMS(c6, strand6, rl6)

		# Expected output
		eqs5, eqe5, eal5 = 42, 418, 377
		eqs6, eqe6, eal6 = 2647, 3725, 1079

		self.assertEqual(qs5, eqs5)
		self.assertEqual(qe5, eqe5)
		self.assertEqual(al5, eal5)
		self.assertEqual(qs6, eqs6)
		self.assertEqual(qe6, eqe6)
		self.assertEqual(al6, eal6)


	def test_cigar2posSMD(self):
		# Test the function cigar2posSMD
		c7, strand7, rl7 = "37727S3378M13D", '-', 41105  
		c8, strand8, rl8 = "1067S7449M11D", '+', 8516  
        
		qs7, qe7, al7 = cigar_parsing.cigar2posSMD(c7, strand7, rl7)
		qs8, qe8, al8 = cigar_parsing.cigar2posSMD(c8, strand8, rl8)

		# Expected output
		eqs7, eqe7, eal7 = 0, 3377, 3391
		eqs8, eqe8, eal8 = 1067, 8515, 7460

		self.assertEqual(qs7, eqs7)
		self.assertEqual(qe7, eqe7)
		self.assertEqual(al7, eal7)
		self.assertEqual(qs8, eqs8)
		self.assertEqual(qe8, eqe8)
		self.assertEqual(al8, eal8)


	def test_cigar2posMDS(self):
		# Test the function cigar2posMDS
		c9, strand9, rl9 = "565M1D623S", '-', 1188   
		c10, strand10, rl10 = "2181M9D4901S", '+', 7082   
        
		qs9, qe9, al9 = cigar_parsing.cigar2posMDS(c9, strand9, rl9)
		qs10, qe10, al10 = cigar_parsing.cigar2posMDS(c10, strand10, rl10)

		# Expected output
		eqs9, eqe9, eal9 = 623, 1187, 566
		eqs10, eqe10, eal10 = 0, 2180, 2190

		self.assertEqual(qs9, eqs9)
		self.assertEqual(qe9, eqe9)
		self.assertEqual(al9, eal9)
		self.assertEqual(qs10, eqs10)
		self.assertEqual(qe10, eqe10)
		self.assertEqual(al10, eal10)


	def test_cigar2posSMDS(self):
		# Test the function cigar2posSMDS
		c11, strand11, rl11 = "2698S1057M28D4430S", '-', 8185  
		c12, strand12, rl12 = "7864S7515M260D38S", '+', 15417   
        
		qs11, qe11, al11 = cigar_parsing.cigar2posSMDS(c11, strand11, rl11)
		qs12, qe12, al12 = cigar_parsing.cigar2posSMDS(c12, strand12, rl12)

		# Expected output
		eqs11, eqe11, eal11 = 4430, 5486, 1085
		eqs12, eqe12, eal12 = 7864, 15378, 7775

		self.assertEqual(qs11, eqs11)
		self.assertEqual(qe11, eqe11)
		self.assertEqual(al11, eal11)
		self.assertEqual(qs12, eqs12)
		self.assertEqual(qe12, eqe12)
		self.assertEqual(al12, eal12)


	def test_cigar2posSMI(self):
		# Test the function cigar2posSMI
		c13, strand13, rl13 = "3282S43752M1269I", '-', 48303   
		c14, strand14, rl14 = "1214S44215M126I", '+', 45555    
        
		qs13, qe13, al13 = cigar_parsing.cigar2posSMI(c13, strand13, rl13)
		qs14, qe14, al14 = cigar_parsing.cigar2posSMI(c14, strand14, rl14)

		# Expected output
		eqs13, eqe13, eal13 = 0, 45020, 43752
		eqs14, eqe14, eal14 = 1214, 45554, 44215

		self.assertEqual(qs13, eqs13)
		self.assertEqual(qe13, eqe13)
		self.assertEqual(al13, eal13)
		self.assertEqual(qs14, eqs14)
		self.assertEqual(qe14, eqe14)
		self.assertEqual(al14, eal14)


	def test_cigar2posMIS(self):
		# Test the function cigar2posMIS
		c15, strand15, rl15 = "3758M19I17947S", '-', 21724 
		c16, strand16, rl16 = "454M3I855S", '+', 1312 
        
		qs15, qe15, al15 = cigar_parsing.cigar2posMIS(c15, strand15, rl15)
		qs16, qe16, al16 = cigar_parsing.cigar2posMIS(c16, strand16, rl16)

		# Expected output
		eqs15, eqe15, eal15 = 17947, 21723, 3758
		eqs16, eqe16, eal16 = 0, 456, 454

		self.assertEqual(qs15, eqs15)
		self.assertEqual(qe15, eqe15)
		self.assertEqual(al15, eal15)
		self.assertEqual(qs16, eqs16)
		self.assertEqual(qe16, eqe16)
		self.assertEqual(al16, eal16)


	def test_cigar2posSMIS(self):
		# Test the function cigar2posSMIS
		c17, strand17, rl17 = "18S528M1I4301S", '-', 4848   
		c18, strand18, rl18 = "36S1199M23I1327S", '+', 2585    
        
		qs17, qe17, al17 = cigar_parsing.cigar2posSMIS(c17, strand17, rl17)
		qs18, qe18, al18 = cigar_parsing.cigar2posSMIS(c18, strand18, rl18)

		# Expected output
		eqs17, eqe17, eal17 = 4301, 4829, 528
		eqs18, eqe18, eal18 = 36, 1257, 1199

		self.assertEqual(qs17, eqs17)
		self.assertEqual(qe17, eqe17)
		self.assertEqual(al17, eal17)
		self.assertEqual(qs18, eqs18)
		self.assertEqual(qe18, eqe18)
		self.assertEqual(al18, eal18)


	def test_alignment_from_satags(self):
		
		# Test the function alignment_from_satags
		sal1 = ['chr1,204201693,-,24154S24767M26I34S,60,589', 'chr1,204201685,+,24831S23987M610D163S,60,3969']
		rl1 = 48981
		sal2 = ['chr20,46787280,+,27S2201M12I1538S,60,87', 'chr20,46788830,+,3028S651M3I96S,45,40', 
			'chr20,46788838,+,2818S642M2D318S,17,46', 'chr14,105563460,+,1963S1051M31I733S,0,358']
		rl2 = 3778

		# Expected output
		eqint1 = [[34, 24826], [24831, 48817]]
		erint1 = [['chr1', 204226458, 204201692, '-'], ['chr1', 204201684, 204226280, '+']]
		equal1 = [60, 60]
		eqint2 = [[27, 2239], [1963, 3044], [2818, 3459], [3028, 3681]]
		erint2 = [['chr20', 46787279, 46789479, '+'], ['chr14', 105563459, 105564509, '+'], 
			['chr20', 46788837, 46789480, '+'], ['chr20', 46788829, 46789479, '+']]
		equal2 = [60, 0, 17, 45]
  
		qint1, rint1, qual1 = cigar_parsing.alignment_from_satags(sal1, rl1)
		qint2, rint2, qual2 = cigar_parsing.alignment_from_satags(sal2, rl2)
		
		self.assertEqual(qint1, eqint1)
		self.assertEqual(rint1, erint1)
		self.assertEqual(qual1, equal1)
		self.assertEqual(qint2, eqint2)
		self.assertEqual(rint2, erint2)
		self.assertEqual(qual2, equal2)


if __name__ == "__main__":
	unittest.main()
