import torch
import numpy as np
from melanoma_image_classifier import Net

def main():
    img_size = 50

    net = Net()
    net.load_state_dict(torch.load('saved_model.pth'))
    net.eval()

    testing_data = np.load("melanoma_testing_data.npy", allow_pickle=True)    

    test_X = torch.Tensor( [item[0] for item in testing_data])
    test_X = test_X / 255

    test_y = torch.Tensor( [item[1] for item in testing_data])

    correct = 0
    total = 0

    with torch.no_grad():
        for i in range(len(test_X)):
            output = net(test_X[i].view(-1, 1, img_size, img_size))[0]
            if output[0]>=output[1]:
                guess = "benign"
            else:
                guess = "melanoma"

            real_label = test_y[i]

            if real_label[0] >= real_label[1]:
                real_class = "benign"
            else:
                real_class = "melanoma"

            if guess == real_class:
                correct +=1
            
            total +=1

    print("Accuracy: ", round(correct/total,3))

if __name__ == "__main__":
    main()