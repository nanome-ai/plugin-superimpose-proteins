import os
import unittest
import sys

plugin_dir = f'{os.getcwd()}/plugin/'
sys.path.append(plugin_dir)

test_directory = 'tests'
suite = unittest.TestLoader().discover(test_directory)

output = unittest.TextTestRunner(verbosity=1).run(suite)

if output.wasSuccessful():
    sys.exit(0)
else:
    sys.exit(1)
