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

qasmGateToProjectQ = {
"u3"   : X,
"u"    : X,
"u2"   : X,
"rx"   : X,
"ry"   : X,
"u1"   : X,
"p"    : X,
"rz"   : X,
"u0"   : X,
"id"   : X,
"t"    : X,
"tdg"  : X,
"s"    : X,
"sdg"  : X,
"z"    : Z,
"x"    : X,
"y"    : Y,
"h"    : H,
"sx"   : X,
"sxdg" : X,
}

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
    qubitIndex   = gQuantumCircuit.find_bit(qubit).registers[0][1]

    print(gQuantumCircuit.find_bit(qubit).registers[0][1])
    print("QUBIT          : ", qubit)
    print("Qubit Reg Name : " + qubitRegName)
    print("Qubit Index    : " + str(qubitIndex))
    print("One Qubit Gate!")

    qasmGateToProjectQ[gateName] | gRegisters[qubitRegName][qubitIndex] 
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

observedString = ""                                                                                     
for regGroup in gRegisters:                                                                             
  print("")                                                                                             
  print("REGISTER GROUP: ", regGroup)                                                                   
  for index in range(0, len(gRegisters[regGroup])):                                                     
    observedString += str(int(gRegisters[regGroup][index]))                                             
    print(f"Index = {index}, value = {int(gRegisters[regGroup][index])}")                            
print("Observed String: ", observedString)                                                           









