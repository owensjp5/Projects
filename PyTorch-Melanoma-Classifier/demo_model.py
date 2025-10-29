import cv2
import numpy as np
import torch
from melanoma_image_classifier import Net
import warnings
warnings.filterwarnings("ignore")

def apply_model(path):
    img_size = 50

    img = cv2.imread(path, cv2.IMREAD_GRAYSCALE)

    img = cv2.resize(img, (img_size,img_size))

    img_array = np.array(img)

    img_array = img_array / 255

    img_array = torch.Tensor(img_array)

    net = Net()
    net.load_state_dict(torch.load("saved_model.pth"))
    net.eval()

    net_out = net(img_array.view(-1, 1, img_size, img_size))[0]

    if net_out[0] >= net_out[1]:
        print()
        print("Prediction: BENIGN")
        print(f"Confidence: {round(float(net_out[0]),3)}")
        print()
    else:
        print()
        print("Prediction: MELANOMA")
        print(f"Confidence: {round(float(net_out[1]),3)}")
        print()