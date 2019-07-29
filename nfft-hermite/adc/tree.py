'''
Defines a basic binary tree class for use in adc.py
'''

class Tree:
	"Simple binary tree."
	
	def __init__(self, data = None, left = None, right = None):
		self.data = data
		self.left = left
		self.right = right