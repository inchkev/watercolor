from PIL import Image, ImageOps
import numpy as np

im = Image.open("paper.jpeg", "r")
im_grayscale = ImageOps.grayscale(im)

array = np.array([[im_grayscale.getpixel((i,j)) for j in range(200)] for i in range(200)], dtype=float)
array -= 50
array /= array.max()
# print(array)

grad = np.gradient(array)
print("dhx = {", end="")
for i in grad[0]:
    print("{",end="")
    for j in i:
        print(f"{j:0.2f},", end="")
    print("},",end="")
print("};")
print("dhy = {", end="")
for i in grad[1]:
    print("{",end="")
    for j in i:
        print(f"{j:0.2f},", end="")
    print("},",end="")
print("};")
