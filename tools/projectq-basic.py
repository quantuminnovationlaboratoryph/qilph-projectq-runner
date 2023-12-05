from projectq          import MainEngine
from projectq.backends import Simulator
from projectq.ops      import (
    CNOT,
    QFT,
    All,
    Barrier,
    BasicMathGate,
    C,
    Entangle,
    H,
    Measure,
    Ph,
    QubitOperator,
    Rx,
    Ry,
    Rz,
    S,
    SqrtSwap,
    SqrtX,
    Swap,
    T,
    Tensor,
    TimeEvolution,
    Toffoli,
    X,
    Y,
    Z,
    get_inverse,
)
 
import sys
from qiskit import *

#=======================================#

gQasmFilename = ""

for parameter in sys.argv:
  if parameter.endswith(".qasm"):
    gQasmFilename = parameter

#=======================================#

gRegisters      = {}
simulator1      = Simulator()
engine1         = MainEngine(backend=simulator1)

gQuantumCircuit = QuantumCircuit.from_qasm_file(gQasmFilename)
gateCounts      = 0
gCircuitQubits  = gQuantumCircuit.num_qubits
gCircuitDepth   = gQuantumCircuit.depth()
gCircuitGates   = gQuantumCircuit.count_ops()




print("")
print("======================================================================")
print("INFO: QUBITS: " + str(gCircuitQubits))
print("INFO: DEPTH : " + str(gCircuitDepth))


print("")
print("======================================================================")
print(gQuantumCircuit)


print("")
print("======================================================================")
print("INFO: QREGS: ", gQuantumCircuit.qregs)
for qreg in gQuantumCircuit.qregs:
  qregSize = qreg.size
  qregName = qreg.name
  gRegisters[qregName] = engine1.allocate_qureg(qregSize)
  print(qreg)



for reg in gRegisters:
  print(reg, gRegisters[reg])



print("")
print("======================================================================")

qasmOneQubitGateLabels = [ 
"u3"  ,
"u"   ,
"u2"  ,
"rx"  ,
"ry"  ,
"u1"  ,
"p"   ,
"rz"  ,
"u0"  ,
"id"  ,
"t"   , 
"tdg" ,
"s"   ,
"sdg" ,
"z"   ,
"x"   ,
"y"   ,
"h"   ,
"sx"  ,
"sxdg",
]

for gate in gQuantumCircuit.data:
  gateName      = gate.operation.name
  gateNumQubits = gate.operation.num_qubits
  gateNumClbits = gate.operation.num_clbits

  #print(gate)
  print("")
  print("Gate Name      : " + gateName)
  print("No. of Qubits  : " + str(gateNumQubits))
  print("No. of Clbits  : " + str(gateNumClbits))

  if gateName.lower() in qasmOneQubitGateLabels:
    qubit        = gate.qubits[0] 
    qubitRegName = gQuantumCircuit.find_bit(qubit)[1][0][0].name
    qubitIndex   = gQuantumCircuit.find_bit(qubit)[0]

    print("QUBIT          : ", qubit)
    print("Qubit Reg Name : " + qubitRegName)
    print("Qubit Index    : " + str(qubitIndex))
    print("One Qubit Gate!")

    #H | gRegisters[qubitRegName][qubitIndex] 
  else:
    print("Unknown 1-Qubit Gate.")

'''
  if gateName.lower() == "cx" or gateName.lower == "cnot":
    controlQubit        = gate.qubits[0] 
    controlQubitRegName = gQuantumCircuit.find_bit(controlQubit)[1][0][0].name
    controlQubitIndex   = gQuantumCircuit.find_bit(controlQubit)[0]

    targetQubit         = gate.qubits[1]
    targetQubitRegName  = gQuantumCircuit.find_bit(targetQubit)[1][0][0].name
    targetQubitIndex    = gQuantumCircuit.find_bit(targetQubit)[0]
    
    print("CX GATE")
    print("CONTROL : ", controlQubitRegName, controlQubitIndex)
    print("TARGET  : ", targetQubitRegName , targetQubitIndex )
    C(X) | (gRegisters[controlQubitRegName][controlQubitIndex], gRegisters[targetQubitRegName][targetQubitIndex])
  elif gateName.lower() == "h":
'''  

print("")
print("======================================================================")

for regGroup in gRegisters:
  All(Measure) | gRegisters[regGroup]
engine1.flush()

print("\n==============================")
observedString = ""
for regGroup in gRegisters:
  print("")
  print("REGISTER GROUP: ", regGroup)
  for index in range(0, len(gRegisters[regGroup])):
    observedString += str(int(gRegisters[regGroup][index]))
    print(f"Index = {index}, value = {int(gRegisters[regGroup][index])}")
   # print("value ", int(gRegisters[regGroup][reg]))
#    print(f"GROUP: {regGroup}, INDEX: {reg}, VALUE: ", gRegisters[regGroup])
#    print(f"GROUP: {regGroup}, INDEX: {reg}, VALUE: ", type(gRegisters[regGroup]))
#    print(f"GROUP: {regGroup}, INDEX: {reg}, VALUE: ", int(gRegisters[regGroup][0]))
#    print(f"GROUP: {regGroup}, INDEX: {reg}, VALUE: ", int(gRegisters[regGroup][1]))
print("Observed String: ", observedString)
print("==============================\n")

#=======================================#

engine0 = MainEngine()
qbit1   = engine0.allocate_qubit()

H       | qbit1
Measure | qbit1

engine0.flush()

print("\n==============================")
print(f"qbit1 type: {type(qbit1)}")
print(f"qbit1  int: {int(qbit1)}")
print("==============================\n")


#=======================================#

engine0 = MainEngine()
qbit1   = engine0.allocate_qubit()
ampl_a0  = engine0.backend.get_amplitude('0', qbit1)
ampl_b0  = engine0.backend.get_amplitude('1', qbit1)
H       | qbit1
ampl_a1  = engine0.backend.get_amplitude('0', qbit1)
ampl_b1  = engine0.backend.get_amplitude('1', qbit1)
engine0.flush()
ampl_a2  = engine0.backend.get_amplitude('0', qbit1)
ampl_b2  = engine0.backend.get_amplitude('1', qbit1)
Measure | qbit1
ampl_a3 = engine0.backend.get_amplitude('0', qbit1)
ampl_b3 = engine0.backend.get_amplitude('1', qbit1)

print("\n==============================")
print(f"qbit1 type: {type(qbit1)}")
print(f"qbit1  int: {int(qbit1)}")
print(f"ampl_a0  : {ampl_a0}")
print(f"ampl_b0  : {ampl_b0}")
print(f"ampl_a1  : {ampl_a1}")
print(f"ampl_b1  : {ampl_b1}")
print(f"ampl_a2  : {ampl_a2}")
print(f"ampl_b2  : {ampl_b2}")
print(f"ampl_a3  : {ampl_a3}")
print(f"ampl_b3  : {ampl_b3}")
print("==============================\n")


#=======================================#

qubitCount = 3
simulator1 = Simulator()
engine1    = MainEngine(backend=simulator1)
qreg       = engine1.allocate_qureg(qubitCount)

All(H)       | qreg
All(Measure) | qreg
engine1.flush()

print("\n==============================")
print(f"qreg     type: {type(qreg)}")
print(f"qreg[0]  type: {type(qreg[0])}")
print(f"qreg[0]   int: {int(qreg[0])}")
print(f"qreg[1]   int: {int(qreg[1])}")
print(f"qreg[2]   int: {int(qreg[2])}")
print("==============================\n")


#============================#
# Engine = Simulator Backend #
#============================#

#eng = MainEngine(backend=UnitarySimulator())

eng = MainEngine(backend=Simulator())
# This is the same as MainEngine()

# With specified randomization seed or repeatability
#eng = MainEngine(backend=Simulator(rnd_seed=10))


simulator1 = Simulator()
simulator2 = Simulator(rnd_seed=10)

engine0    = MainEngine()
engine1    = MainEngine(backend=simulator1)
engine2    = MainEngine(backend=simulator2)

#========#
# Qubits #
#========#

# qubit = eng.allocate_qubit()  # allocate 1 qubit
qubitCount = 5
qreg       = engine2.allocate_qureg(qubitCount)

print(qreg)



#======================#
# Operator Application #
#======================#

# H | qubit  # apply a Hadamard gate
# NOT | qubit
# Measure | qubit  # measure the qubit

H | qreg[0]
H | qreg[1]
H | qreg[2]
H | qreg[3]
H | qreg[4]

#All(H)       | qreg
All(Measure) | qreg


eng.flush()  # flush all gates (and execute measurements)
#print(f"Measured {int(qubit)}")  # output measurement result
print("QREG: ", int(qreg[0]))
print("QREG: ", int(qreg[1]))
print("QREG: ", int(qreg[2]))
print("QREG: ", int(qreg[3]))
print("QREG: ", int(qreg[4]))
print("QREG TYPE: ", type(qreg[4]))
print(f"Measured {int(qreg[0])}")  # output measurement result
print(f"Measured {int(qreg[1])}")  # output measurement result
print(f"Measured {int(qreg[2])}")  # output measurement result
print(f"Measured {int(qreg[3])}")  # output measurement result
print(f"Measured {int(qreg[4])}")  # output measurement result



