#!/usr/bin/env python

import unittest

class TestPrimers(unittest.TestCase):

    def setUp(self):
        # add primers to database
        raise NotImplementedError

    def test_retrival(self):
        # get
        raise NotImplementedError

    def test_design(self):
        # get with design
        raise NotImplementedError

    def test_isupper(self):
        self.assertTrue('FOO'.isupper())
        self.assertFalse('Foo'.isupper())

    def test_split(self):
        s = 'hello world'
        self.assertEqual(s.split(), ['hello', 'world'])
        # check that s.split fails when the separator is not a string
        with self.assertRaises(TypeError):
            s.split(2)

if __name__ == '__main__':
    unittest.main()
