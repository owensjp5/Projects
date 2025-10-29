import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from melanoma_image_classifier import Net

def main():
    img_size = 50

    training_data = np.load("melanoma_training_data.npy", allow_pickle=True)

    train_X = torch.Tensor( [item[0] for item in training_data]  )
    train_X = train_X / 255

    train_y = torch.Tensor( [item[1] for item in training_data]  )

    net = Net()

    optimizer = optim.Adam(net.parameters(), lr=0.001)

    loss_func = nn.MSELoss()

    batch_size = 100

    epochs = 10

    for epoch in range(epochs):

        for i in range(0,len(train_X), batch_size):

            print(f"EPOCH {epoch+1}, % complete: {i/len(train_X)}")

            batch_X = train_X[i: i+batch_size].view(-1, 1, img_size, img_size)
            batch_y = train_y[i: i+batch_size]

            optimizer.zero_grad()

            outputs = net(batch_X)

            loss = loss_func(outputs,batch_y)

            loss.backward()

            optimizer.step()

    torch.save(net.state_dict(), "saved_model.pth")

if __name__ == "__main__":
    main()