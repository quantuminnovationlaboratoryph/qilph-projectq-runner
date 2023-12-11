import math
import numpy as np
from projectq          import MainEngine
from projectq.backends import Simulator
from projectq.ops      import *
'''
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
'''
 
import sys
from qiskit import *

#=======================================#

debug = False

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

qasmMultiQubitGateLabels = [ 
"cx"  ,
"cu"  ,
]

# This function takes a small threshold number and a complex number cNum = p + qi, where a and b are 
# the two real components of cNum, and return new complex number cNum' = r + si where r is equal to
# p if the absolute value of p's fractional component is greater than the threshold value otherwise
# r is the integral part of p and s is equal to q if the absolute value of q's fractional component
# is greater than the threshold value otherwise s is the integral part of q.

def zeroDownFrac (cNum, threshold):
  realComp = np.modf(np.real(cNum))
  realInt  = realComp[1]
  realFrac = realComp[0]
  imagComp = np.modf(np.imag(cNum))
  imagInt  = imagComp[1]
  imagFrac = imagComp[0]
  realFrac = 0.0 if np.fabs(realFrac) < threshold else realFrac
  imagFrac = 0.0 if np.fabs(imagFrac) < threshold else imagFrac
  return complex(realInt + realFrac, imagInt + imagFrac)

def qasmU3Matrix ( _t, _p, _l, thres):
  _k = _p + _l
  p  = np.cos(_t/2.0)
  q  = np.sin(_t/2.0)
  e0 = 1
  e1 = -(np.cos(_l) + np.sin(_l)*0j)
  e2 =  (np.cos(_p) + np.sin(_p)*0j)
  e3 =  (np.cos(_k) + np.sin(_k)*0j)
  w  = p * e0
  x  = q * e1
  y  = q * e2
  z  = p * e3
  ww = zeroDownFrac(w, thres)
  xx = zeroDownFrac(x, thres)
  yy = zeroDownFrac(y, thres)
  zz = zeroDownFrac(z, thres)
  if debug == True:
    print("theta   =", _t)
    print("phi     =", _p)
    print("lambda  =", _l)
    print("p (cos) = ", p)
    print("q (sin) = ", q)
    print("e0      = ", e0)
    print("e1      = ", e1)
    print("e2      = ", e2)
    print("e3      = ", e3)
    print("w       = ", w)
    print("x       = ", x)
    print("y       = ", y)
    print("z       = ", z)
    print("ww      = ", ww)
    print("xx      = ", xx)
    print("yy      = ", yy)
    print("zz      = ", zz)
  return np.array([[ww, xx],[yy, zz]])


thres = 0.0000000001
#print(H.matrix)
#print(myGate.matrix)
print("Hadmard Matrix from U3:\n", qasmU3Matrix(math.pi/2.0, 0, math.pi, thres), "\n")
print("NOT Matrix from U3:\n", qasmU3Matrix(math.pi, 0, math.pi, thres), "\n")

idGate = MatrixGate([[1,0],
                     [0,1]])

SqrtXdag = get_inverse(SqrtX)

qasmGateToProjectQ = {
"u3"      : MatrixGate,
"u"       : MatrixGate,
"u2"      : MatrixGate,
"rx"      :         Rx,
"ry"      :         Ry,
"u1"      : MatrixGate,
"p"       :         Ph,
"rz"      :         Rz,
"u0"      :     idGate, #
"id"      :     idGate, #
"t"       :          T, 
"tdg"     :       Tdag, 
"s"       :          S, 
"sdg"     :       Sdag,
"z"       :          Z, 
"x"       :          X, 
"y"       :          Y, 
"h"       :          H, 
"sx"      :      SqrtX,
"sxdg"    :   SqrtXdag, #
"barrier" :    Barrier,
"measure" :    Measure,
}

gate = MatrixGate([[0, 1], [1, 0]])

for gate in gQuantumCircuit.data:
  gateName      = gate.operation.name.lower()
  gateParams    = gate.operation.params
  gateNumQubits = gate.operation.num_qubits
  gateNumClbits = gate.operation.num_clbits


  # For Debugging
  if gateName in ["barrier", "measure"]:
    continue

  print("\n")
  print("Gate Name      : " + gateName)
  print("Gate Parameters: " + str(gateParams))
  print("No. of Qubits  : " + str(gateNumQubits))
  print("No. of Clbits  : " + str(gateNumClbits))
  
  if gateName in qasmOneQubitGateLabels:

    qubit        = gate.qubits[0] 
    qubitRegName = gQuantumCircuit.find_bit(qubit)[1][0][0].name
    qubitIndex   = gQuantumCircuit.find_bit(qubit).registers[0][1]

    print("Qubit Reg Name : " + qubitRegName)
    print("Qubit Index    : " + str(qubitIndex))
    
    if gateName.lower() in ["u3", "u", "u2", "u1"]:
      if   len(gateParams) == 3:
         gateParams2 = (gateParams[0], gateParams[1], gateParams[2], thres)
      elif len(gateParams) == 2:
         gateParams2 = (  math.pi/2.0, gateParams[0], gateParams[1], thres)
      else:
         gateParams2 = ( 0000000000.0, 00000000000.0, gateParams[0], thres)
      matrixParam = qasmU3Matrix(*gateParams2)
      qasmGateToProjectQ[gateName](matrixParam)   | gRegisters[qubitRegName][qubitIndex] 
    elif gateName in ["rx", "ry", "rz", "p"]:
      qasmGateToProjectQ[gateName](gateParams[0]) | gRegisters[qubitRegName][qubitIndex] 
    else:
      qasmGateToProjectQ[gateName]                | gRegisters[qubitRegName][qubitIndex] 

    print("One Qubit Gate!")

  else:
    continue
    #print("Unknown 1-Qubit Gate.")

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









