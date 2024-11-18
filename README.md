Project: Segmentation of abdominal aorta when there is a contrast agent.

Overview

Main goal of this project was the segmentation and the 3D reconstruction of the abdominal aorta which is extracted from CT data.

The method is divided into two parts. In the first part we estimate the position and the dimension of the aortic lumen using state-of-the-art object tracking techniques. The second part employs curve fitting methods in order to detect the boundaries of the aortic lumen with accuracy.

In particular, the proposed method uses the Kalman Filter to track the aortic cross-section in consecutive CT images. The observations needed by  the Kalman procedure are extracted with the Circle Hough Transformation, based on the assumption that the morphological structure of the aortic cross-section is approximately a circle. A robust Level Set method is then applied to compensate the approximation error and efficiently estimate the cross-section.

The algorithms and the mathematical tools developed during the project prove feasibility for an accurate and reliable method for the segmentation of the abdominal aneurysm from CT data, that in the future could be used to benefit patients with aortic aneurysms.

Method

The ﬁrst part is based on a “rough” estimation of the position of the aortic lumen.  For this purpose we assume that the aortic-cross sections are circular and we try to track the position of the aorta using the Kalman Filter. The observation (measurement) needed by the Kalman procedure is extracted from the Hough Circle Transformation (HCT). In the ﬁrst image we do not have any previous knowledge about the cross-section position so we apply the HCT in the extended region of the image. This way we extract the initial state vector (one of the initial conditions) needed for the Kalman Filtering. Then the ﬁlter estimates and corrects HCT measurements and moreover deﬁnes smaller regions (Region of Interest-ROI) for the HCT application (gif. ﬁgure Kalman Tracking ). When the recursive procedure ends a level set method is applied.
