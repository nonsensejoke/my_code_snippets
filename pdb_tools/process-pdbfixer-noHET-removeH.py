from pdbfixer import PDBFixer
from openmm.app import PDBFile, Modeller
from openmm.app.element import hydrogen


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

modeller = Modeller(fixer.topology, fixer.positions)
hydrogens = [atom for atom in modeller.topology.atoms() if atom.element.symbol == "H"]
modeller.delete(hydrogens)

# 注意：我们在这里故意不调用 fixer.addMissingHydrogens()，因此不会添加氢原子 [8, 9]

# 验证删除后是否还有 H

remaining_h = []

for atom in modeller.topology.atoms():

    if atom.element == hydrogen or atom.name.strip().upper().startswith("H"):

        remaining_h.append(atom)

print("Remaining possible H:", len(remaining_h))

# 5. 保存结果文件
PDBFile.writeFile(modeller.topology, modeller.positions, open(sys.argv[2], 'w'), keepIds=True) # [3, 10]
