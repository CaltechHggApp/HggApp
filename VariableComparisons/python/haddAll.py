from samples import sampleMap
import os

stub = """
hadd -f /tmp/%s.root output/%s/*.root; mv /tmp/%s.root output/
"""

for k,v in sampleMap.iteritems():
    print k
    cmd = stub % (k,k,k,)
    print cmd
    os.system(cmd)
