#!/usr/bin/python

# dot -Tsvg -o out.svg in.gv
# digraph cdeprov2dot {
# graph [rankdir = "RL"];
# node [fontname="Helvetica" fontsize="8" style="filled" margin="0.0,0.0"];
# edge [fontname="Helvetica" fontsize="8"];
# "nnnn1" [label="nnnn1" URL="url1" shape="box" fillcolor="lightsteelblue1"]
# "nnnn2" [label="nnnn2" URL="url2" shape="box" fillcolor="lightsteelblue1"]
# "nnnn1" -> "nnnn2" [label="" color="blue"]
# }

from __future__ import division
import commands
import re
import os, sys, time, datetime, glob, json
import argparse
from anytree import Node, RenderTree

#####################################  Need a Tree? ########################################

class Library():
	
	## Python
	python = Node("python")
	
	Keras = Node('Keras', parent = python)
	pandas = Node("pandas", parent=python)
	
	abc = Node("abc", parent = python)									
	asyncio = Node("asyncio", parent = python)		
	backcall = Node("backcall", parent = python)									
	beautifulsoup4 = Node("beautifulsoup4", parent = python)									
	blinker = Node("blinker", parent = python)									
	certifi = Node("certifi", parent = python)									
	cftime = Node("cftime", parent = python)									
	chardet = Node("chardet", parent = python)									
	collections = Node('collections', parent = python)									
	concurrent = Node("concurrent", parent = python)		
	cpython = Node("cpython", parent = python)									
	ctypes = Node('ctypes', parent = python)									
	dateutil = Node("dateutil", parent=pandas)									
	distutils = Node('disutils', parent = python)									
	email = Node('email', parent = python)									
	encodings = Node("encodings", parent = python)									
	h5py = Node('h5py', parent = Keras)									
	hs_restclient = Node('hs_restclient', parent = python)									
	html = Node("html", parent = python)		
	http = Node("http", parent = python)									
	idna = Node('idna', parent = python)									
	importlib = Node("importlib", parent = python)									
	ipykernel = Node('ipykernel', parent = python)									
	ipyleaflet = Node('ipyleaflet', parent = python)									
	IPython = Node('IPython', parent = python)									
	ipython_genutils = Node("ipython_genutils", parent = python)									
	ipywidgets = Node('ipywidgets', parent = python)									
	jedi = 	Node('jedi', parent = python)								
	jinja2 = Node("jinja2", parent = python)									
	json = Node("json", parent = python)		
	jupyter_client = Node("jupyter_client", parent = python)							
	jupyter_core = Node("jupyter_core", parent = python)						
	logging = Node("logging", parent = python)	
	markupsafe = Node("markupsafe", parent = python)									
	matplotlib = Node("matplotlib", parent = python)	
	multiprocessing = Node("multiprocessing", parent = python)									
	netCDF4 = Node("netCDF4", parent = python)									
	numpy = Node('numpy', parent = h5py)									
	numpy = Node('numpy', parent = Keras)									
	numpy = Node("numpy", parent=pandas)									
	numpy = Node("numpy", parent=python)									
	oauthlib = Node("oauthlib", parent = python)																		
	parso = Node("parso", parent = python)									
	pexpect = Node("pexpect", parent = python)									
	PIL = Node("PIL", parent = python)									
	pkg_resources = Node("pkg_resources", parent = python)			### Different on PyPI						
	prompt_toolkit = Node("prompt_toolkit", parent = python)									
	protobuf = Node("protobuf", parent = python)									
	ptyprocess = Node("ptyprocess", parent = python)	
	pyasn1 = Node("pyasn1", parent = python)	
	pycurl = Node("pycurl", parent = python)	
	pydoc_data = Node("pydoc_data", parent = python)								
	pygments = Node("pygments", parent = python)									
	pygobject =	Node("pygobject", parent = python)
	pyparsing = Node("pyparsing", parent = python)									
	pyproj = Node("pyproj", parent = python)	
	pysumma = Node("pysumma", parent = python)		### Not on Pypi
	pytz = Node("pytz", parent=pandas)									
	requests = Node("requests", parent = python)									
	requests_toolbelt = Node("requests-toolbelt", parent=python) 								
	six = Node('six', parent = h5py)									
	sqlite3	= Node('sqlite3', parent = python)	
	tkinter = Node("tkinter", parent = python)									
	tornado = Node("tornado", parent = python)									
	traitlets = Node("traitlets", parent = python)									
	traittypes = Node("traittypes", parent = python)									
	unittest = Node("unittest", parent = python)																	
	urllib3 = Node("urllib3", parent = python)									
	wcwidth = Node("wcwidth", parent = python)									
	xarray = Node("xarray", parent = python)	
	xml = Node("xml", parent = python)								
	zmq = Node("zmq", parent = python)									

	
	## R
	R = Node('R')

	A3 = Node('A3', parent = R)
	abbyyR = Node('abbyyR', parent = R)
	abc = Node('abc', parent = R)	
	assertthat = Node('assertthat', parent = R)	
	bitops = Node('bitops', parent = R)
	caTools = Node('caTools', parent = R)
	colorspace = Node('colorspace', parent = R)
	curl = Node('curl', parent = R)	
	digest = Node('digest', parent = R)
	gdata = Node('gdata', parent = R)
	data__table = Node('data.table', parent = R)
	ggplot2 = Node('ggplot2', parent = R)
	MASS = Node('MASS', parent = R)
	randomForest = Node('randomForest', parent = R)
	ROCR = Node('ROCR', parent = R)
	RSocrata = Node('RSocrata', parent = R)
	gplots = Node('gplots', parent = R)
	gtable = Node('gtable', parent = R)
	gtools = Node('gtools', parent = R)
	httr = Node('httr', parent = R)
	jsonlite = Node('jsonlite', parent = R)
	KernSmooth = Node('Kernsmooth', parent = R)
	labeling = Node('labeling', parent = R)
	lazyeval = Node('lazyeval', parent = R)
	mime = Node('mime', parent = R)
	munsell = Node('munsell', parent = R)
	plyr = Node('plyr', parent = R)
	R6 = Node('R6', parent = R)
	Rcpp = Node('Rcpp', parent = R)	
	RUnit = Node('RUnit', parent = R)	
	scales = Node('scales', parent = R)
	tibble = Node('tibble', parent = R)
	
#for pre, fill, node in RenderTree(Library.python):
#	print("%s%s" % (pre, node.name))

#for pre, fill, node in RenderTree(Library.R):
#	print("%s%s" % (pre, node.name))

#print(python3.children)

############################################################################################

## Jason: Storing Variables For Documentation Write-Up
inPut = []			# To print inputs (as shown on provenance graph)
processes = []			# To print processes (as shown on provenance graph)
outPut = []			# To print outputs (as shown on provenance graph)
dependencies = []		# To print all filtered system dependencies
unused_dependencies = []
piddict = {}			# Match execute PID and source
reqpackages = set()		# Pick out required package from libaries list
pkgcheck = set()
unwantedlist = ['cpython','cache','/dev/','/var/', '/etc/', '/sys/','/proc/', '/usr/share/','/tmp/', '/sudo/', 'linux','/locale/']
unwantedpaths = set()
pkginfo = []

libdict = {key.replace('__','.'):None for key, value in Library.__dict__.items() if not key.startswith('__') and not callable(key)}
libdict.pop('python', None)
libdict.pop('R', None)


#Filter Package for PIP SEARCH
def func1(s):
    temp = s.split('lib/')
    if len(temp)==1:
        return False
    temp = temp[1].split('/',1)
    if len(temp)==1:
        return False
    else:    
        return False if '/' in temp[1] else temp[1]
    
def func2(s):
    temp = s.split('site-packages/')
    if len(temp)==1:
        return False
    return False if '/' in temp[1] else temp[1]
    
def func3(s):
    temp = s.split('dist-packages/')
    if len(temp)==1:
        return False
    return False if '/' in temp[1] else temp[1]

func_set = [func1, func2, func3]

def test_path(func_set, path):
    for func in func_set:
        res = func(path)
        if res:
            return res
    else:
        return False


def isFilteredPath(path):
  if re.match('\/proc\/', path) is None \
    and re.match('.*\/lib\/', path) is None \
    and re.match('\/etc\/', path) is None \
    and re.match('\/var\/', path) is None \
    and re.match('\/dev\/', path) is None \
    and re.match('\/sys\/', path) is None \
    and re.match('.*\/R\/x86_64-pc-linux-gnu-library\/', path) is None \
    and re.match('.*\/usr\/share\/', path) is None:
    return False
  else: 
    return True
  
parser = argparse.ArgumentParser(description='Process provenance log file.')
parser.add_argument('--nosub', action="store_true", default=False)
parser.add_argument('--nofilter', action="store_true", default=False)
parser.add_argument('--withfork', action="store_true", default=False)
parser.add_argument('-f', action="store", dest="fin_name", default="provenance.cde-root.1.log")
parser.add_argument('-d', action="store", dest="dir_name", default="./gv")
parser.add_argument('--withgraph', action="store_true", default=False)

args = parser.parse_args()

showsub = not args.nosub
filter = not args.nofilter
withfork = args.withfork
dir = args.dir_name
logfile = args.fin_name
withgraph = args.withgraph
#Hai: commented the line 49
#pickle_message_generator.deliver_messages(logfile) #delivers messages about missing pickle files

meta = {}
colors=['yellow','pink','lightgreen','lightblue'] 
colorid=0
# BIG TODO: colors are currently arranged to make the later overcolors the former
#    so that sub-graph overcolors the parent graph nodes
re_set = {}

def processMetaData(line):
  m = re.match('# @(\w+): (.*)$', line)
  (key, value) = m.group(1, 2)
  meta[key] = value
  if key == "namespace":
    re_set['rmns'] = re.compile('"[^"]*' + meta['namespace'] + '\/')
  return key

def makePathNode(path):
  filename=os.path.basename(path).replace('"','\\"')
  node=path.replace('\\', '\\\\').replace('"','\\"')
  nodedef='"' + node + '"[label="' + filename + '", shape="", fillcolor=' + colors[colorid] + ', tooltip="' + node + '"]'
  return (node, nodedef)
  
def getEdgeStr(node, pathnode, deptype):
  if deptype == "wasGeneratedBy":
    return '"' + pathnode + '" -> ' + node
  elif deptype == "used" or deptype == "rw":
    return node + ' -> "' + pathnode + '"'
  else:
    return 'ERROR -> ERROR'

def printArtifactDep(node, path, deptype, filter):





  #Jason: Record filtered path as system dependencies
	if (filter or isFilteredPath(path)):		# Mark path that uses libraries. Otherwise assign system to path
		flag = 0		
		for word in unwantedlist:
			if word in path:
				unwantedpaths.add(path)
				flag = 1

		if flag == 0:
			pkgcheck.add(test_path(func_set, path))			
			for key in libdict.keys():
				if (('site-packages/' + key) in path) or (('library/' + key) in path) or (('dist-packages/' + key) in path) or ('/usr/lib/' in path and '/'+key in path):
					if 'dist-info' in path or 'egg-info' in path or 'nspkg' in path:
						unused_dependencies.append(piddict.get(words[1]) + ' --> ' + path)
				 	elif (key + ' --> ' + path) not in dependencies:						
						dependencies.append(key + ' --> ' + path)
						reqpackages.add(key)
					flag = 1
		
		if flag == 0:
			pkgcheck.add(test_path(func_set, path))
			if ('site-packages/py' in path) or ('library/py' in path) or ('dist-packages/py' in path):		
				if 'dist-info' in path or 'egg-info' in path or 'nspkg' in path:
						unused_dependencies.append(piddict.get(words[1]) + ' --> ' + path)		
				elif ('py --> ' + path) not in unused_dependencies:
					dependencies.append('py --> ' + path)
					reqpackages.add('py')
				flag = 1

		if flag == 0:
			pkgcheck.add(test_path(func_set, path))
			if 'urllib' in path:
				if 'dist-info' in path or 'egg-info' in path or 'nspkg' in path:
						unused_dependencies.append(piddict.get(words[1]) + ' --> ' + path)		
				elif ('urllib --> ' + path) not in unused_dependencies:
					dependencies.append('urllib --> ' + path)
					reqpackages.add('urllib')
				flag = 1

		if flag == 0:
			pkgcheck.add(test_path(func_set, path))			
			try:	
				if (piddict.get(words[1]) + ' --> ' + path) not in unused_dependencies:
					unused_dependencies.append(piddict.get(words[1]) + ' --> ' + path)
			except:
				pass





	if (not filter or not isFilteredPath(path)):
		direction = 'dir="both" ' if deptype == "rw" else ''
		(pnode, pdef) = makePathNode(path)
		edge = getEdgeStr(node, pnode, deptype)
		fout.write(pdef + '\n')
		fout.write(edge + ' [' + direction +'label="' + deptype + '" color=""];\n') # blue
		if showsub:
			try:# TOFIX: what about spawn node
				pid_graph[node].write(pdef + '\n')
				pid_graph[node].write(edge + ' [' + direction +'label="' + deptype + '" color=""];\n') # blue
			except:
				pass
		if withgraph:
			fgraph.write(re_set['rmns'].sub('"/', edge) + '\n')



## Jason: Record Inputs and Outputs
		pnodename = pnode.split('/')
		if deptype == 'used':
			if (path) not in inPut:
				inPut.append(path)

		if deptype == 'wasGeneratedBy':
			if (path) not in outPut:
				outPut.append(path)

		if deptype == 'rw':
			if (path) not in inPut:
				inPut.append(path)
			if (path) not in outPut:
				outPut.append(path)



# open input output files
try:
  fin = open(logfile, 'r')
except IOError:
  print "Error: can\'t find file " + logfile + " or read data\n"
  sys.exit(-1)

while 1:
  line = fin.readline()
  if re.match('^#.*$', line) or re.match('^$', line):
    if re.match('^# @.*$', line):
      processMetaData(line)
    continue
  else:
    break

# prepare graphic directory
if not os.path.exists(dir):
  os.makedirs(dir)
os.system("rm -f " + dir + "/*.gnu " + dir + "/*.svg " + dir + "/*.gv " + dir + "/*.html")

fout = open(dir + '/main.gv', 'w')
fout.write("""digraph G {
graph [rankdir = "RL" ];
node [fontname="Helvetica" fontsize="8" style="filled" margin="0.0,0.0"];
edge [fontname="Helvetica" fontsize="8"];
"cdenet" [label="XXX" shape="box" fillcolor=""" +colors[colorid]+ "];\n" + \
"\"namespace" + meta['fullns'] + '"[shape=box label="' + meta['agent'] + "@" + meta['fullns'] + '" color=' +colors[colorid]+ ']')
f2out = open(dir + '/main.process.gv', 'w')
f2out.write("""digraph cdeprovshort2dot {
graph [rankdir = "RL" ];
node [fontname="Helvetica" fontsize="8" style="filled" margin="0.0,0.0"];
edge [fontname="Helvetica" fontsize="8"];
"cdenet" [label="XXX" shape="box" fillcolor=""];
"unknown" [label="unknown" shape="box" fillcolor=""];
""") # blue lightsteelblue1
fhtml = open(dir + '/main.html', 'w')
fhtml.write("""<h1>Overview</h1>
<a href=main.svg>Full Graph</a><br/>
<a href=main.process.svg>Process Graph</a><br/>
""")
fhtml.close()

if withgraph:
  fgraph = open(dir + '/main.graph', 'w')

active_pid = {'0':'unknown'}
info_pid = {}
cde_pid = -1
pid_desc = {'cdenet':'[label="XXX" shape="box" fillcolor="" URL="main.process.svg"]', # blue
  'unknown':'[label="unknown" shape="box" fillcolor="" URL="main.process.svg"]'} # blue
iodpid_desc = {'cdenet':'[label="XXX" shape="box" fillcolor="" URL="main.process.svg"]', # blue
  'unknown':'[label="unknown" shape="box" fillcolor="" URL="main.process.svg"]'}
pid_starttime = {}
pid_mem = {}
pid_graph = {}
counter = 1

while 1:
  if re.match('^#.*$', line) or re.match('^$', line):
    if re.match('^# @.*$', line):
      if processMetaData(line) == 'fullns':
        colorid += 1
        fout.write('"namespace' + meta['fullns'] + '"[shape=box label="' + meta['agent'] + "@" + meta['fullns'] + '" color=' +colors[colorid]+ ']\n')
    line = fin.readline()
    continue
  line = line.rstrip('\n').replace('\\', '\\\\').replace('"','\\"')
  words = line.split(' ', 6)
  pid = words[1]
  action = words[2]
  path = '' if len(words) < 4 else words[3]
  path = path.replace('"', '\"')
  
  if pid in active_pid:
    node = active_pid[pid] # possible NULL
  else:
    cde_pid = pid
    active_pid[pid] = 'cdenet'
  
  if action == 'EXECVE': # this case only, node is the child words[3]
    nodename = words[3] + '_' + str(counter)
    node = '"' + nodename + '"'
    label = 'PID: ' + words[3] + "\\n" + os.path.basename(words[4])
    iodlabel = os.path.basename(words[4])
    iodpath = words[4]
    piddict.update({words[3]:os.path.basename(words[4])})
    #print(piddict)

    title = time.ctime(int(words[0])) + ' ' + words[4] + ' ' + ''.join(words[5:])
    #    .replace(', \\"',', \\n\\"').replace('[\\"','\\n[\\"')
    counter += 1
    active_pid[words[3]] = node # store the dict from pid to unique node name
    info_pid[words[3]] = {'path': words[4], 'dir': words[5]}
    
    # main graph
    pid_desc[node] = ' [label="' + label + '" tooltip="' + title + '" shape="box" fillcolor="' + colors[colorid] + '" URL="' + nodename + '.prov.svg"]' # lightsteelblue1
    iodpid_desc[node] = iodlabel
    fout.write(node + pid_desc[node] + '\n')
    fout.write(node + ' -> ' + active_pid[pid] + ' [label="wasTriggeredBy" color=""]\n') # darkblue

    if withgraph:
      fgraph.write(node + ' : ' + json.dumps(info_pid[words[3]]) + '\n')
      fgraph.write(node + ' -> ' + active_pid[pid] + '\n')
    
    # main process graph
    f2out.write(node + pid_desc[node] + '\n')
    f2out.write(node + ' -> ' + active_pid[pid] + ' [label="wasTriggeredBy" color=""]\n') # darkblue
    
    # html file of this pid
#     htmlf = open('gv/' + nodename + '.html', 'w')
#     htmlf.write('Memory Footprint: <img src=' + nodename + '.mem.svg /><br/>' +
#       'Provenance Graph: <object data=' + nodename + '.prov.svg type="image/svg+xml"></object>')
#     htmlf.close()

    if showsub:
      # mem graph of this pid
      pid_mem[node] = open(dir + '/' + nodename + '.mem.gnu', 'w')
      pid_mem[node].write('set terminal svg\nset output "' + nodename + '.mem.svg"\n' +
        'set xlabel "Time (seconds) from ' + time.ctime(int(words[0])) + '\n' +
        'set ylabel "Memory Usage (kB)"\n' +
        'set xrange [-1:]\n' +
        'plot "-" title "' + label + '" with lines\n' +
        '0\t0\n')
      pid_starttime[node] = int(words[0])
      
      parentnode = active_pid[pid]
      # prov graph of this pid
      pid_graph[node] = open(dir + '/' + nodename + '.prov.gv', 'w')
      pid_graph[node].write('digraph "' + nodename + '" {\n'+
        """graph [rankdir = "RL" ];
        node [fontname="Helvetica" fontsize="8" style="filled" margin="0.0,0.0"];
        edge [fontname="Helvetica" fontsize="8"];
        """ +
        'Memory [label="" shape=box image="' + nodename + '.mem.svg"];\n' +
        node + ' [label="' + label + '" shape="box" fillcolor=""]\n' + #green
        parentnode + pid_desc[parentnode] + '\n' +
        node + ' -> ' + parentnode + ' [label="wasTriggeredBy" color=""]\n') #darkblue



## Jason: Records Processes (except shell files connected to 'XXX', marked as input)
      if parentnode == 'cdenet' and '.sh' in iodpath:	
        inPut.append(iodpath)
      if (iodpath) in inPut:
				pass
      else:
				if iodpath not in processes:
					processes.append(iodpath)



      # prov graph of parent pid
      try:
        pid_graph[parentnode].write(
          node + pid_desc[node] + '\n' +
          node + ' -> ' + parentnode + ' [label="wasTriggeredBy" color=""]\n') #darkblue
      except:
        pass

  elif action == 'SPAWN':
    node = '"' + words[3] + '_' + str(counter) + '"'
    label = 'fork ' + words[3]
    counter += 1
    parentnode = active_pid[pid]
    active_pid[words[3]]=node # store the dict from pid to unique node name
    info_pid[words[3]]=info_pid[words[1]]
    pid_desc[node] = pid_desc[parentnode]
    
    if withfork: # not handle the withgraph case
      # main graph
      fout.write(node + ' [label="' + label + '" shape="box" fillcolor=""];\n') #azure
      fout.write(node + ' -> ' + parentnode + ' [label="wasTriggeredBy" color=""];\n') #darkblue
      
      # main process graph
      f2out.write(node + ' [label="' + label + '" shape="box" fillcolor=""];\n') #azure
      f2out.write(node + ' -> ' + parentnode + ' [label="wasTriggeredBy" color=""];\n') #darkblue
    else:
      active_pid[words[3]]=active_pid[pid]

  elif action == 'EXIT':
    if showsub:
      try:
        pid_mem[node].close()
        del pid_mem[node]
        pid_graph[node].write("}")
        pid_graph[node].close()
        del pid_graph[node]
      except:
        pass
      del active_pid[pid]
    
  elif action == 'READ':
    pinfo = info_pid[words[1]]
    if os.path.abspath(pinfo['path']) != path:
      printArtifactDep(node, path, "used", filter)
    
  elif action == 'WRITE':
    printArtifactDep(node, path, "wasGeneratedBy", filter)
  elif action == 'READ-WRITE':
    printArtifactDep(node, path, "rw", filter)
    
  elif action == 'MEM':
    if showsub:
      rel_time = int(words[0]) - pid_starttime[node]
      try:
        pid_mem[node].write(str(rel_time) + '\t' + str(int(words[3])>>10) + '\n')
      except KeyError:
        pass
  
  line = fin.readline()
  if line == '':
    break
  
fout.write("}")
fout.close()
if withgraph:
  fgraph.close()
f2out.write("}")
f2out.close()
fin.close()

if showsub:
  for node in pid_graph:
    try:
      pid_graph[node].write("}")
      pid_graph[node].close()
    except:
      pass

  for node in pid_mem:
    try:
      pid_mem[node].close()
    except:
      pass

def removeMultiEdge(filename):
  lines = [line for line in open(filename)]
  newlines = lines[:4]
  newlines = newlines + list(set(lines[4:-1])) + ['}']
  f = open(filename, 'w')
  for line in newlines:
    f.write(line)
  f.close()

# covert created graphviz and gnuplot files into svg files
os.chdir(dir)
os.system("dot -Tsvg main.process.gv -o main.process.svg")
print("Processing main.gv, this might take a long time ...\n")
removeMultiEdge("main.gv")
os.system("dot -Tsvg main.gv -o main.svg")
#os.system("unflatten -f -l3 main.gv | dot -Tsvg -Edir=none -Gsplines=ortho -o main.svg")
os.system("gnuplot *.gnu")
for fname in glob.glob("./*.prov.gv"):
  removeMultiEdge(fname)
  os.system("dot -Tsvg " + fname + " -o " + fname.replace('prov.gv', 'prov.svg'))
os.system('echo "<h1>Processes</h1>" >> main.html')
os.system('ls *.prov.svg | while read l; do echo "<a href=$l>$l</a><br/>" >> main.html; done')
os.system('echo "<h1>Memory Footprints</h1>" >> main.html')
os.system('ls *.mem.svg | while read l; do echo "<a href=$l>$l</a><br/>" >> main.html; done')


# Handling PIP SEARCH Packages
nspkg = set()
extras = []

try:	
	pkgcheck.remove(False)
	pkgcheck.remove('lib-dynload')
	pkgcheck.remove('dist-packages')
	pkgcheck.remove('site-packages')
except:
	pass


for i in pkgcheck:
	if 'info' in i:
		pkginfo.append(i)
	if 'nspkg' in i:
		nspkg.add(i)
	if '.py' in i:
		extras.append(i.split('.')[0])

for i in nspkg:
	pkgcheck.remove(i)
for i in pkginfo:
	pkgcheck.remove(i)
for i in extras:
	pkgcheck.add(i)




pkginfo = [i.replace(i, i.replace('.dist-info','').replace('.egg-info','')) for i in pkginfo]

def comparver(ver1,ver2):
	comp1 = ver1.split('.')
	comp2 = ver2.split('.')
	for i in range(min(len(comp1),len(comp2))):
		if int(comp1[i]) > int(comp2[i]):
			return ver1
		elif int(comp1[i]) < int(comp2[i]):
			return ver2
	if len(comp1) > len(comp2):
		return ver1
	else:
		return ver2

pkginfo = sorted(pkginfo)
newpkginfo = []

if len(pkginfo) == 0:
	pass
else:
	prevword=(pkginfo[0].split('-'))[0]
	prevver=(pkginfo[0].split('-'))[1]

	for i in range(1, len(pkginfo)):
		word = (pkginfo[i].split('-'))[0]
		ver = (pkginfo[i].split('-'))[1]
		if word != prevword:
			newpkginfo.append(prevword+'{'+prevver+'}')
			prevword = word
			prevver = ver
		elif word == prevword:
			prevword = word
			prevver = comparver(prevver, ver)

	newpkginfo.append(prevword+"{"+prevver+"}")

pkgtoremove = set()
pkgtoappend = set()
for i in pkgcheck:
	for j in newpkginfo:
		if i.lower() in j.lower():
			pkgtoremove.add(i)
			pkgtoappend.add(j)


for i in pkgtoremove:
	pkgcheck.remove(i)
for i in pkgtoappend:
	pkgcheck.add(i)




hyphenlst = []
identifiedpkg = []
for i in list(pkgcheck):
	if '_' in i:
		hyphenlst.append(i.replace('_','-'))
		pkgcheck.remove(i)

for i in list(pkgcheck):		
	if commands.getoutput('pip search '+i):			
		identifiedpkg.append(i)

for i in hyphenlst:
	if commands.getoutput('pip search '+i):			
		identifiedpkg.append(i)






## Jason: Documenting IOD
writeIOD = open('IOD.txt', 'w')
writeIOD.write('Documentation of Project belonging to: ' + meta['agent'] + "@" + meta['fullns'])
writeIOD.write('\n' + '-'*60 + '\n\n')

writeIOD.write('INPUTS')
writeIOD.write('\n' + '-'*6 + '\n')
for i in sorted(inPut):
	writeIOD.write(i + '\n')
	
writeIOD.write('\nPROCESSES')
writeIOD.write('\n' + '-'*9 + '\n')
for i in sorted(processes):
	writeIOD.write(i + '\n')

writeIOD.write('\nOUTPUTS')
writeIOD.write('\n' + '-'*7 + '\n')
for i in sorted(outPut):
	writeIOD.write(i + '\n')

writeIOD.write('\nREQUIRED PACKAGES')
writeIOD.write('\n' + '-'*17 + '\n')
writeIOD.write('No. of packages required: ' + str(len(reqpackages)) + '\n')
writeIOD.write('Packages: \n')
for i in sorted(reqpackages):
  writeIOD.write(i + ' ')
writeIOD.write('\n\nNo. of packages required through PIP SEARCH: ' + str(len(identifiedpkg)) + '\n')
writeIOD.write('Packages: \n')
for i in sorted(identifiedpkg):
  writeIOD.write(i + ' ')
#writeIOD.write('\n\nAccuracy: '+ str(round((len(reqpackages)/len(identifiedpkg))*100, 2)) + '%')

#writeIOD.write('\n\nPACKAGE INFO')
#writeIOD.write('\n' + '-'*12 + '\n')
#for i in sorted(newpkginfo):
#	writeIOD.write(i + '\n')

writeIOD.write('\n\nDEPENDENCIES')
writeIOD.write('\n' + '-'*12 + '\n')
for i in sorted(dependencies):
	writeIOD.write(i + '\n')

writeIOD.write('\nUNUSED DEPENDENCIES')
writeIOD.write('\n' + '-'*20 + '\n')
for i in sorted(unused_dependencies):
	writeIOD.write(i + '\n')

writeIOD.write('\nUNWANTED PATH')
writeIOD.write('\n' + '-'*20 + '\n')
for i in sorted(unwantedpaths):
	writeIOD.write(i + '\n')

writeIOD.close()
