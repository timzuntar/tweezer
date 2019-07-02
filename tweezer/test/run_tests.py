import unittest

diagnostics = unittest.TestLoader().discover(".")

unittest.TextTestRunner().run(diagnostics)