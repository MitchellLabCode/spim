import json
"""
Specialize to 2 views (0, 90 deg)
576 lines tall at 250 fps --> 1ms exp, 3ms delay time
576 lines in width or so, as needed, check 90 deg view
500 reps
30 min
1-2 micron step dz

"""


def create_array(zmin, zmax, dz):
    return [zmin + i * dz for i in range(int((zmax - zmin) / dz) + 1)]


with open('E:/jbutler/test500fps/add_stack_fromJan_101524.json', 'r') as file:
  contents = json.load(file)

dz = 2  # microns, the step in z for each view
minz = [-85, -302]  # the starting value of z for each view (STACK 0, stack 1)
maxz = [230, 0]   # the final value of z for each view... Note: ranges are in columns. ex, stack0 z range -154,-422
rot_angles = [0, 270]  # the first angle rotation
# rot_incrememnt = 90  # can be set to whatever value you want
# rot_angles = [rot_angle, rot_angle + rot_incrememnt]
xvals = [4874, 4874]
yvals = [-150, -254]

stacks = contents['data']['stacks']

dup_stack = stacks[-1]
zvals = [create_array(minz[0], maxz[0]+dz, dz), create_array(minz[1], maxz[1]+dz,  dz)]

# for each view
for j in [0, 1]:
    rot_angle = rot_angles[j]
    # add a new "stack" with fixed z
    for i in range(0, len(zvals[j])):
        print(zvals[j][i])
        new_stack = json.loads(json.dumps(dup_stack))
        new_stack['name'] = f"stack_{i + j * len(zvals[0])}"
        z_value = zvals[j][i]
        for element in new_stack['elements']:
            if element['name'] == 'z':
                element['start'] = z_value
                element['end'] = z_value
            if element['name'] == 'r':
                element['start'] = rot_angle
                element['end'] = rot_angle
            if element['name'] == 'x':
                element['start'] = xvals[j]
                element['end'] = xvals[j]
            if element['name'] == 'y':
                element['start'] = yvals[j]
                element['end'] = yvals[j]

        stacks.append(new_stack)


with open('/sid/add_multipleStacks_RotateAngle_20250227_v3.json', 'w') as f_out:
    json.dump(contents, f_out, indent=2)

print('done')
print(contents)