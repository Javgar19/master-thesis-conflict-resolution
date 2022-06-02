from sacred import Experiment

ex = Experiment("hello_cs")

@ex.config
def cfg(_log):
	message = "This is printed by function {}"

@ex.capture
def foo(message):
	print(message.format("foo"))

@ex.automain
def main(message):
	foo()
	foo("Overriding message")