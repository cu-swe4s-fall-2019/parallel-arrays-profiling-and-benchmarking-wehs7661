import unittest
import plot_gtex


class TestSearchingMethods(unittest.TestCase):
    def test_linear_search(self):
        """ This function tests the functionality of the linear search method
        """
        L = [1, 2, 3, 4, 5]
        r = plot_gtex.linear_search(3, L)
        self.assertEqual(r, 2)

        r = plot_gtex.linear_search(100, L)
        self.assertEqual(r, -1)

    def test_linear_search_hits(self):
        """ This function tests the functionality of the linear search method
        """
        L = ['a', 'b', 'b', 'a', 'a', 'a', 'b', 'a', 'a']

        r = plot_gtex.linear_search_all_hits('a', L)
        self.assertEqual(r, [0, 3, 4, 5, 7, 8])

        r = plot_gtex.linear_search_all_hits('c', L)
        self.assertEqual(r, [])

    def test_binary_search(self):
        """ This function tests the functionality of the binary search method
        """
        L = [[1, 0], [2, 0], [3, 0], [4, 0], [5, 0]]

        for i in range(1, 6):
            r = plot_gtex.binary_search(i, L)
            self.assertEqual(r, 0)

        r = plot_gtex.binary_search(100, L)
        self.assertEqual(r, -1)


if __name__ == '__main__':
    unittest.main()
