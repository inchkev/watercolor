from PIL import Image, ImageOps
import numpy as np

im = Image.open("paper_2.jpeg", "r")
im_grayscale = ImageOps.grayscale(im)

array = np.array([[im_grayscale.getpixel((i*2,j*2)) for j in range(512)] for i in range(512)], dtype=float)
array += 200
array /= array.max()
array *= 254
print(np.min(array))
print(np.max(array))

im = Image.fromarray(array)
im = im.convert('RGB')
im.save("h_2.ppm")
