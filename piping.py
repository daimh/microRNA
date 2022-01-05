import threading, subprocess, os, tempfile, random, sys
class PipingException(Exception): pass

class FuncBase:
	def work(self, line, fout): pass
	def close(self, fout): pass

class Pipe:
	def __init__(self, step, fin=subprocess.DEVNULL):
		self.stdout = fin
		self._steps = []
		self._threads = []
		self._fifos = []
		self.append(step)
		self.stdin = self._steps[0].stdin
	def append(self, step):
		if isinstance(step, list): #step is command
			self._steps.append(subprocess.Popen(step, stdin=self.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True))
			if not isinstance(self.stdout, int):
				self.stdout.close()
			self.stdout = self._steps[-1].stdout
			thr = threading.Thread(target=self._get_stderr_thread, args=(self._steps[-1],))
			self._threads.append(thr)
			thr.start()
		elif isinstance(step, FuncBase): #step is Function
			self._steps.append(step)
			fifo_out = self._mkfifo()
			if isinstance(self.stdout, int):
				if self.stdout != subprocess.PIPE:
					raise PipingException('FuncBase pipe must have subprocess.PIPE or a file like object as input')
				fin = self._mkfifo()
				stdin_is_fifo = True
			else:
				fin = self.stdout
				stdin_is_fifo = False
			thr = threading.Thread(target=self._append_func_thread, args=(step, fin, fifo_out))
			self._threads.append(thr)
			thr.start()
			self.stdout = open(fifo_out)
			if stdin_is_fifo:
				step.stdin = open(fin, 'w')
			else:
				step.stdin = fin
		else:
			raise PipingException('Only pipe.FuncBase or a list of string can be used as parameter')
	def _append_func_thread(self, func, fin, fifo_out):
		fout = open(fifo_out, 'w')
		if isinstance(fin, str):
			fin = open(fin)
		for ln in fin: func.work(ln, fout)
		fin.close()
		func.close(fout)
		fout.close()
	def _get_stderr_thread(self, step):
			step.serr = step.stderr.read()
			step.stderr.close()
	def _mkfifo(self):
		tmp = tempfile.gettempdir() + '/piping.' + ''.join([random.choice('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ') for i in range(16)])
		if os.path.exists(tmp): raise PipingException(f'A temporary file with the same name {tmp} already exists')
		os.mkfifo(tmp)
		self._fifos.append(tmp)
		return tmp
	def tee(self, pipes):
		thread = threading.Thread(target=self._tee_thread, args=(pipes, ))
		thread.start()
		self._threads.append(thread)
	def _tee_thread(self, pipes):
		for p in pipes:
			if p.stdin is None: raise PipingException('One pipe in tee is not initliazed with subprocess.PIPE')
		for ln in self.stdout:
			for p in pipes:
				print(ln, end='', file=p.stdin)
		self.stdout.close()
		for p in pipes:
			p.stdin.close()
	def close(self):
		for thr in self._threads: thr.join()
		for step in self._steps:
			if not isinstance(step, subprocess.Popen): continue
			step.wait()
			if step.returncode != 0:
				print('#', step.args, 'returned', step.returncode, 'with following stderr\n', step.serr, file=sys.stderr)
		self.stdout.close()
		for tmp in self._fifos: os.unlink(tmp)
