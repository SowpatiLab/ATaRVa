import ATARVA.pyssw as pyssw
import ATARVA.ssw_lib as ssw_lib
import ctypes as ct
import os, sys


def suppress_stderr():
    """Redirects stderr to /dev/null to suppress C++ error messages."""
    sys.stderr.flush()  # Flush any buffered stderr data
    stderr_fileno = sys.stderr.fileno()  # Get stderr file descriptor
    devnull = os.open(os.devnull, os.O_WRONLY)  # Open /dev/null for writing
    os.dup2(devnull, stderr_fileno)  # Redirect stderr to /dev/null
    os.close(devnull)  # Close descriptor

def restore_stderr(original_stderr):
    """Restores the original stderr after suppression."""
    sys.stderr.flush()
    os.dup2(original_stderr, sys.stderr.fileno())
    os.close(original_stderr)
    # sys.stderr = sys.__stderr__  # Restore original stderr

class Inputs:
    def __init__(self, target, query):
        # self.sLibPath = 'Complete_Striped_Smith_Waterman_Library/src'  # Set your libssw.so path
        self.sLibPath = False # Set your libssw.so path
        self.nMatch = 2
        self.nMismatch = 2
        self.nOpen = 3
        self.nExt = 1
        self.bProtein = False
        self.sMatrix = ''  # Set your matrix file if needed
        self.bPath = True  # False
        self.nThr = 0
        self.bBest = False
        self.bSam = True  # False
        self.bHeader = False
        self.target = target  # Set your target file path
        self.query = query  # Set your query file path

def stripSW(args):
    lEle = []
    dEle2Int = {}
    dInt2Ele = {}
    
    lEle = ['A', 'C', 'G', 'T', 'N']
    for i,ele in enumerate(lEle):
        dEle2Int[ele] = i
        dEle2Int[ele.lower()] = i
        dInt2Ele[i] = ele
    nEleNum = len(lEle)
    lScore = [0 for i in range(nEleNum**2)]
    for i in range(nEleNum-1):
        for j in range(nEleNum-1):
            if lEle[i] == lEle[j]:
                lScore[i*nEleNum+j] = args.nMatch
            else:
                lScore[i*nEleNum+j] = -args.nMismatch
    
    ssw = ssw_lib.CSsw(args.sLibPath)
    
    sRSeq = args.target
    sQSeq = args.query
    
    qNum = pyssw.to_int(sQSeq, lEle, dEle2Int)
    mat = (len(lScore) * ct.c_int8) ()
    mat[:] = lScore
    qProfile = ssw.ssw_init(qNum, ct.c_int32(len(sQSeq)), mat, len(lEle), 2)
    nMaskLen = len(sQSeq) // 2
    nFlag = 2
    
    rNum = pyssw.to_int(sRSeq, lEle, dEle2Int)
    if nMaskLen<15:
        original_stderr = os.dup(sys.stderr.fileno())
        suppress_stderr()
        res = pyssw.align_one(ssw, qProfile, rNum, len(sRSeq), args.nOpen, args.nExt, nFlag, nMaskLen)
        restore_stderr(original_stderr)
    else:
        res = pyssw.align_one(ssw, qProfile, rNum, len(sRSeq), args.nOpen, args.nExt, nFlag, nMaskLen)
    sCigar, sQ, sA, sR = pyssw.buildPath(sQSeq, sRSeq, res[4], res[2], res[8])
    return sA, [res[4], res[5]]
