from pdbfixer import PDBFixer
from openmm.app import PDBFile

import sys

# 1. 加载要处理的文件
fixer = PDBFixer(filename=sys.argv[1]) # [3]

# remove certain chain if necessary 
#fixer.removeChains(chainIds=['A'])

# 2. 识别缺失的残基
fixer.findMissingResidues() # [5]
# 可选操作：如果您不想补全缺失的整段残基（只补全现有残基缺失的原子），可以加上下面这行代码清空记录
fixer.missingResidues = {} # [6]

fixer.removeHeterogens(keepWater=False)

# 3. 识别缺失的重原子
fixer.findMissingAtoms() # [7]

# 4. 实际添加所有缺失的重原子
fixer.addMissingAtoms() # [8]

# 注意：我们在这里故意不调用 fixer.addMissingHydrogens()，因此不会添加氢原子 [8, 9]

# 5. 保存结果文件
PDBFile.writeFile(fixer.topology, fixer.positions, open(sys.argv[2], 'w'), keepIds=True) # [3, 10]
