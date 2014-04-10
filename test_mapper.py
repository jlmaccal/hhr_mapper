#!/usr/bin/env python
# encoding: utf-8


import unittest
import mapper


class TestFileIO(unittest.TestCase):
    def setUp(self):
        self.lines = open('sequence.hhr').readlines()
        self.line_ranges = mapper.find_line_ranges(self.lines)

    def test_find_line_ranges_gives_correct_ranges(self):
        self.assertEqual(self.line_ranges[0], (114, 140))
        self.assertEqual(self.line_ranges[1], (140, 166))

    def test_find_line_ranges_gives_correct_number_of_matches(self):
        self.assertEqual(len(self.line_ranges), 104)

    def test_parse_match_number_is_correct(self):
        match_number = mapper.parse_match_number(
            self.lines[self.line_ranges[0][0]])
        self.assertEqual(match_number, 1)

    def test_parse_description_gives_correct_information(self):
        pdb, chain, desc = mapper.parse_description(
            self.lines[self.line_ranges[0][0] + 1])
        self.assertEqual(pdb, '4fr9')
        self.assertEqual(chain, 'A')
        self.assertIn('Uncharacterized', desc)

    def test_gather_query_sequence_lines_gets_correct_lines(self):
        lines = mapper.get_query_sequence_lines(
            self.lines[self.line_ranges[0][0]:self.line_ranges[0][1]])
        self.assertEqual(len(lines), 2)
        self.assertTrue(lines[0].startswith('Q T0644'))
        self.assertTrue(lines[1].startswith('Q T0644'))

    def test_gather_template_sequence_lines_gets_correct_lines(self):
        lines = mapper.get_template_sequence_lines(
            self.lines[self.line_ranges[0][0]:self.line_ranges[0][1]])
        self.assertEqual(len(lines), 2)
        self.assertTrue(lines[0].startswith('T 4fr9_A'))
        self.assertTrue(lines[1].startswith('T 4fr9_A'))

    def test_extract_residue_numbers_has_correct_start_and_end(self):
        lines = mapper.get_query_sequence_lines(
            self.lines[self.line_ranges[0][0]:self.line_ranges[0][1]])
        start1, end1 = mapper.extract_residue_range_from_sequence_line(lines[0])
        start2, end2 = mapper.extract_residue_range_from_sequence_line(lines[1])
        self.assertEqual(start1, 30)
        self.assertEqual(end1, 109)
        self.assertEqual(start2, 110)
        self.assertEqual(end2, 166)

    def test_extract_sequence_from_sequence_line_has_correct_sequence(self):
        lines = mapper.get_query_sequence_lines(
            self.lines[self.line_ranges[0][0]:self.line_ranges[0][1]])
        seq = mapper.extract_sequence_from_sequence_line(lines[0])
        self.assertTrue(seq.startswith('GYL'))
        self.assertTrue(seq.endswith('SYN'))

    def test_gather_query_ss_pred_gets_correct_lines(self):
        lines = mapper.get_query_ss_pred_lines(
            self.lines[self.line_ranges[0][0]:self.line_ranges[0][1]])
        self.assertEqual(len(lines), 2)
        self.assertTrue(lines[0].startswith('Q ss_pred'))
        self.assertTrue(lines[1].startswith('Q ss_pred'))

    def test_gather_tempate_ss_pred_gets_correct_lines(self):
        lines = mapper.get_template_ss_pred_lines(
            self.lines[self.line_ranges[0][0]:self.line_ranges[0][1]])
        self.assertEqual(len(lines), 2)
        self.assertTrue(lines[0].startswith('T ss_pred'))
        self.assertTrue(lines[1].startswith('T ss_pred'))

    def test_gather_tempate_dssp_gets_correct_lines(self):
        lines = mapper.get_template_dssp_lines(
            self.lines[self.line_ranges[0][0]:self.line_ranges[0][1]])
        self.assertEqual(len(lines), 2)
        self.assertTrue(lines[0].startswith('T ss_dssp'))
        self.assertTrue(lines[1].startswith('T ss_dssp'))

    def test_extract_ss_from_ss_line_has_correct_sequence(self):
        lines = mapper.get_query_ss_pred_lines(
            self.lines[self.line_ranges[0][0]:self.line_ranges[0][1]])
        ss = mapper.extract_ss_from_ss_line(lines[0])
        self.assertTrue(ss.startswith('CCC'))
        self.assertTrue(ss.endswith('CCC'))

    def test_get_query_ss_pred_gets_correct_ss_string(self):
        ss = mapper.get_query_ss(
            self.lines[self.line_ranges[0][0]:self.line_ranges[0][1]])
        self.assertEqual(len(ss), 137)
        self.assertTrue(ss.startswith('CCCCCH'))
        self.assertTrue(ss.endswith('CCCCC'))

    def test_get_template_ss_pred_gets_correct_ss_string(self):
        ss = mapper.get_template_ss(
            self.lines[self.line_ranges[0][0]:self.line_ranges[0][1]])
        self.assertEqual(len(ss), 137)
        self.assertTrue(ss.startswith('CccCCH'))
        self.assertTrue(ss.endswith('ccCCCC'))

    def test_get_template_dssp_gets_correct_ss_string(self):
        ss = mapper.get_template_dssp(
            self.lines[self.line_ranges[0][0]:self.line_ranges[0][1]])
        self.assertEqual(len(ss), 137)
        self.assertTrue(ss.startswith('TCCCC'))
        self.assertTrue(ss.endswith('CCCCC'))

    def test_get_query_sequence_is_correct(self):
        seq, start, end = mapper.get_query_sequence(
            self.lines[self.line_ranges[0][0]:self.line_ranges[0][1]])
        self.assertEqual(len(seq), 137)
        self.assertEqual(start, 30)
        self.assertEqual(end, 166)
        self.assertTrue(seq.startswith('GYL'))
        self.assertTrue(seq.endswith('PRV'))

    def test_get_template_sequence_is_correct(self):
        seq, start, end = mapper.get_template_sequence(
            self.lines[self.line_ranges[0][0]:self.line_ranges[0][1]])
        self.assertEqual(len(seq), 137)
        self.assertEqual(start, 8)
        self.assertEqual(end, 144)
        self.assertTrue(seq.startswith('GYL'))
        self.assertTrue(seq.endswith('PRV'))


if __name__ == '__main__':
    unittest.main()
