# BIOT 6900 · Module 10

Sean Cox

In this assignment, I fine-tuned a pretrained DenseNet121 on a binary PathMNIST task: colorectal tumour versus normal tissue. I evaluated the clean model, inspected its predictions with Grad-CAM and occlusion, then tested whether a synthetic red corner tag could become a shortcut. The tag is added to tumour tiles in both the training and validation sets, while tests measure it's effect.

In the saved run, the clean transfer model reached 95.5% validation accuracy and 93.3% test accuracy. The red-tag model reached 99.8% validation accuracy, but tumour recall fell from 100.0% with the tag to 57.7% without it. When the tag was added to normal test tiles, 75.7% were called tumour by the confounded model, compared with 10.7% by the clean control. Grad-CAM and occlusion also assigned more evidence to the tag than its 3.1% image area would suggest. Together, these results support that the model learned the synthetic artifact, while the remaining recall without it shows the model did not rely on the tag exclusively. The experiment does not test robustness across independent scanners, stains, or sites.

I did use AI to help me write part of the code. 