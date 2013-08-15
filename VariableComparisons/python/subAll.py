from samples import sampleMap
import os

stub = """
batchSub.py --submit --DirName=%s --OutputDir=/dev/null --OutputDir=output/%s --FilesPerJob=2 --CMD="sed -e 's:##INPUT##:#IFL#:g' -e 's:##OUTPUT##:#OF1#:g' < run.C > #OF0#" --outputName=run.C --CMD="root -b -q #OF0#" %s
"""

for k,v in sampleMap.iteritems():
    print k
    cmd = stub % (k,k,v)
    print cmd
    os.system(cmd)
