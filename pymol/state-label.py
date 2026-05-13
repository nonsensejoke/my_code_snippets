from pymol import cmd

def make_state_label(sel="aa", x=None, y=None, z=None, obj=""):
    """
    用法：
      make_state_label
      make_state_label aa
      make_state_label aa, protein
      make_state_label aa, protein, 5
      make_state_label none, 20, 15, 5

    说明：
      默认从 selection 'aa' 提取坐标；
      如果不想用 selection，写 none 后面跟 x,y,z。
    """

    # 自动找 trajectory object
    if not obj:
        objs = [o for o in cmd.get_names("objects") if o != "state_label"]
        if not objs:
            print("ERROR: no molecular object found")
            return
        obj = objs[0]

    nstates = cmd.count_states(obj)
    if nstates < 1:
        print(f"ERROR: object/selection '{obj}' has no states")
        return

    # 情况 1：从 selection 提取坐标
    if sel.lower() not in ("none", "manual", "xyz"):
        if cmd.count_atoms(f"({sel})") != 1:
            print(f"ERROR: selection '{sel}' must contain exactly 1 atom")
            print(f"       Current atom count: {cmd.count_atoms(f'({sel})')}")
            return

        xyz = cmd.get_atom_coords(f"({sel})", state=cmd.get_state())
        x, y, z = xyz[0], xyz[1], xyz[2]

    # 情况 2：手动指定 xyz
    else:
        if x is None or y is None or z is None:
            print("ERROR: manual xyz mode needs x, y, z")
            print("Example: make_state_label none, 20, 15, 5")
            return
        x, y, z = float(x), float(y), float(z)

    cmd.delete("state_label")

    for s in range(1, nstates + 1):
        cmd.pseudoatom(
            "state_label",
            pos=[x, y, z],
            state=s,
            label=f"State: {s}"
        )

    cmd.hide("everything", "state_label")
    cmd.show("labels", "state_label")
    cmd.set("label_color", "black", "state_label")
    cmd.set("label_size", 24, "state_label")

    # 如果不想覆盖已有 movie，就注释掉这一行
    cmd.mset(f"1 -{nstates}")

    print(f"state_label created at [{x:.3f}, {y:.3f}, {z:.3f}] using object '{obj}'")

cmd.extend("make_state_label", make_state_label)