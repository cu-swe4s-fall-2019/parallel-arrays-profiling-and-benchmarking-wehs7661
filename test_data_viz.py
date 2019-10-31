import unittest
import data_viz
import os
import random


class TestBoxPlot(unittest.TestCase):
    def test_box_plot(self):
        """This function test if the plot is saved properly
        """
        output_name = 'boxplot_test.png'
        data, data_list = [], []
        meta = ['1', '2', '3', '4', '5']
        for i in range(5):
            for j in range(100):
                data.append(random.randint(1, 1000))
            data_list.append(data)
        check_bf = os.path.exists(output_name)
        data_viz.boxplot(data_list, meta, 'Distribution', 'Value',
                         'Unit test for boxplot function', output_name)
        check_af = os.path.exists(output_name)
        self.assertFalse(check_bf)
        self.assertTrue(check_af)
        os.remove(output_name)


if __name__ == '__main__':
    unittest.main()
