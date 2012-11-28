{
 "metadata": {
  "name": "Dan Domogala (Student Microbiome Project)"
 }, 
 "name": "Dan Domogala (Student Microbiome Project)", 
 "nbformat": 2, 
 "worksheets": [
  {
   "cells": [
    {
     "cell_type": "code", 
     "collapsed": true, 
     "input": "mapping_file_url = \"https://dl.dropbox.com/u/2868868/student_microbiome_project_mapping_31Oct2012.tsv\"\nmapping_file_name = \"student_microbiome_project_mapping_31Oct2012.tsv\"\npersonal_id_field = 3\noutput_fp = '/mnt/ddd38/ipynb_out/abx_disturb_column.txt'", 
     "language": "python", 
     "outputs": [], 
     "prompt_number": 19
    }, 
    {
     "cell_type": "code", 
     "collapsed": false, 
     "input": "!wget $mapping_file_url", 
     "language": "python", 
     "outputs": [
      {
       "output_type": "stream", 
       "stream": "stdout", 
       "text": "--2012-10-31 21:19:34--  https://dl.dropbox.com/u/2868868/student_microbiome_project_mapping_31Oct2012.tsv\nResolving dl.dropbox.com... "
      }, 
      {
       "output_type": "stream", 
       "stream": "stdout", 
       "text": "107.20.198.68, 107.22.246.178, 23.21.176.62, ...\nConnecting to dl.dropbox.com|107.20.198.68|:443... connected."
      }, 
      {
       "output_type": "stream", 
       "stream": "stdout", 
       "text": "HTTP request sent, awaiting response... "
      }, 
      {
       "output_type": "stream", 
       "stream": "stdout", 
       "text": "200 OK\nLength: 2032388 (1.9M) [text/tab-separated-values]\nSaving to: &#96;student_microbiome_project_mapping_31Oct2012.tsv&apos;\n\n\n 0% [                                       ] 0           --.-K/s              "
      }, 
      {
       "output_type": "stream", 
       "stream": "stdout", 
       "text": "\n99% [=====================================&gt; ] 2,031,322   9.64M/s              \n100%[======================================&gt;] 2,032,388   9.63M/s   in 0.2s    \n\n2012-10-31 21:19:34 (9.63 MB/s) - &#96;student_microbiome_project_mapping_31Oct2012.tsv&apos; saved [2032388/2032388]\n"
      }
     ], 
     "prompt_number": 5
    }, 
    {
     "cell_type": "code", 
     "collapsed": true, 
     "input": "antibiotic_disturbed_individuals = [32, 40, 157, 121, 147, 146, 234, 236, 115, 123, 157, 10, 05, 121, 233, 61, 49, 114, 103, 207, 15, 44, 33, 117, 103, 114, 15, 136, 237, 115]", 
     "language": "python", 
     "outputs": [], 
     "prompt_number": 4
    }, 
    {
     "cell_type": "code", 
     "collapsed": false, 
     "input": "output_f = open(output_fp,'w')\nfor line in open(mapping_file_name,'U'):\n    # remove any leading or trailing white space\n    line = line.strip()\n    if not line or line.startswith('#'):\n        # do nothing on blank or header/comment lines\n        pass\n    else:\n        fields = line.split('\\t')\n        current_personal_id = fields[personal_id_field]\n        try:\n            current_personal_id_number = int(current_personal_id[-3:])\n        except ValueError:\n            output_f.write('No\\n')\n        else:\n            if current_personal_id_number in antibiotic_disturbed_individuals:\n                output_f.write('Yes\\n')\n            else:\n                output_f.write('No\\n')\noutput_f.close()\n        ", 
     "language": "python", 
     "outputs": [], 
     "prompt_number": 21
    }, 
    {
     "cell_type": "code", 
     "collapsed": true, 
     "input": "", 
     "language": "python", 
     "outputs": [], 
     "prompt_number": "&nbsp;"
    }
   ]
  }
 ]
}