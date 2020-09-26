from __future__ import print_function

import json
import os.path

class Checkpoint():
	'Class for checkpointing progress'

	def __init__(self, fn):
		self.fn = fn

	def loadJson(self):
		done=dict()
		if os.path.isfile(self.fn):
			print("\n**********************************")
			print("Resuming from previous checkpoint.")
			print("**********************************\n")
			with open(self.fn) as f:
				done = json.load(f)
		return done
	
	def writeJson(self, d):
		with open(self.fn, 'w' ) as json_file:
			json.dump(d, json_file)

				
