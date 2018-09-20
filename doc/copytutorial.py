from __future__ import print_function, unicode_literals
import os
import glob
import shutil

for oldpath in glob.iglob(os.path.join('tutorial', '*.ipynb')):
    oldname = os.path.basename(oldpath)
    newname = oldname.replace(' ', '_')
    newpath = os.path.join('source', 'tutorial', newname)
    shutil.copy(oldpath, newpath)