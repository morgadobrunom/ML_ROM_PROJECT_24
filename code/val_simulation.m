net = importTensorFlowNetwork("model");
net.initialize()
input.BCs = dlarray(Reduced_Data{1},"C")
input.state = dlarray()

F = predict(net,input.BCs, input.state)