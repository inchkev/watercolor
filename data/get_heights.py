from PIL import Image, ImageOps
import numpy as np

im = Image.open("paper_1.jpeg", "r")
im_grayscale = ImageOps.grayscale(im)

array = np.array([[im_grayscale.getpixel((i,j)) for j in range(512)] for i in range(512)], dtype=float)
array += 50
array /= array.max()
array *= 254
print(np.min(array))
print(np.max(array))

im = Image.fromarray(array)
im = im.convert('RGB')
im.save("h_1.ppm")

grad = np.gradient(array)
im = Image.fromarray(grad[0]*2550)
im = im.convert('RGB')
im.save("dhx_1.ppm")

im = Image.fromarray(grad[1]*2550)
im = im.convert('RGB')
im.save("dhy_1.ppm")
# print(array)

# grad = np.gradient(array)
# print("float h[200][200] = {", end="")
# for i in array:
#     print("{",end="")
#     for j in i:
#         print(f"{j:0.4f}".rstrip('0')+",", end="")
#     print("},",end="")
# print("};")
# print("float dhx[200][200] = {", end="")
# for i in grad[0]:
#     print("{",end="")
#     for j in i:
#         print(f"{j:0.4f}".rstrip('0')+",", end="")
#     print("},",end="")
# print("};")
# print("float dhy[200][200] = {", end="")
# for i in grad[1]:
#     print("{",end="")
#     for j in i:
#         print(f"{j:0.4f}".rstrip('0')+",", end="")
#     print("},",end="")
# print("};")
