N1: 
Line 128-129: 
It looks like it's trying too hard. Try something more natural and just say what it is instead of what it is not. Something like: "In addition, to show that there is no significant radiation induced artifact, an indepdent dose-rate test was performed ..."

N2:
Line 135-139: 
I don't understand what you want to say. Please just treat the dilute sample as the same type of sample in the manuscript and stop making a big deal explaining what it is not and what it does not prove. We disclosed the sample concentration, we'll let the readers make their own judgement whether this test suffice or not. 

As a side note: from my observation at 12-ID-B, the state that the damage is most likely to occure is the states in the middle, as the ones shown in the SI. The final arrested state does not have radiation damage. The reason is that g2 is nearly flat. One way to tell if there is radiation damage is to analyze the frames in different sub-sections and see if they overlap. And you can't get a flat g2 unless each portion of it is flat. We can reanalyze the final arrested states by running XPCS correlation over the first 50,000 frames and the second 50,000 frames for all the measurements being averaged, but that takes too much time, and I would rather not do it unless the reviewer asks. So for now just focus on what the data shows and stop talking about what the data does not show.

N3:
Line 141: 
Why do you mention quadrature sum and uncertainty at all? Was there anywhere in main and SI where the uncertainty of SAXS from 12-ID was actually used? If not, please drop it. Stop making the manuscript and SI more precise that necessary.

N4:
Line 189:
Shouldn't it be in the same paragraph as the previous one?

N5:
Line 208: 
Calling it 'Speckle contrast' is confusing, as it is the speckle contrast from the static reference. Rename it something like 'Calibration of beamline coherencen coefficient". 

N6:
Line 209: Refer to Eq. S4 when introducing beta.

N7:
Line 215-225: Again, focus on what the data shows instead of what it does not show. For example line 218: "better-counted bins depart measurably from at and so report residual structure in the standard rather than the instrument." can probably be omitted. Same with line 221, discusssion of contrast outlier measured at a different Q. The contrast is established by the fit in Fig. S4b. Just describe what the data points are, how the fits are done, and why the fit has to cut off higher Q because of longitudinal coherence. The end.

N8:
Line 226:
The check of whether different beta changes the fitting result is unnecessary. As long as the code uses the beta determined in Fig. S4b, just say what beta was used and skip this entire paragraph.

N9:
Line 231:
Similar to N6, refer to Eq. S4. Open the subsection with something more natural.


N10:
Line 256:
Section 5 should be just about structural reversibility and should only contain 5.1 and 5.2. 5.3 is part of the 12-ID data analysis and should belong in 3.1. 5.4 is part of the data set that is shown in 4.1 and 4.3 so they should belong together. In other words, the sections should be:

1. Recombinant Synthesis and Purification of (VPAVG)30. This part can stay the same.
2. Sample Environments and X-ray Measurements. 2.1 and 2.2 can stay. Change the title to something like "Sample Environments and X-ray Measurement Setup".
3. Put 2.3, 'Radiation-damage mitigation and controls', to a new Section 3 because it is important.
4. Section 4 will be the old 3.1 and 5.3.  
5. Section 5 will be the old 3.2, 4.1, 4.2, 4.3, 5.4 and 5.5. 
6. Section 6 will be the old 5.1 and 5.2. The name will be something like Macro and Microscale Structure Reversibility. The section should be dedicated to reversibility.
7. Section 7 and 8 will be the old Section 6 and 7, respectively.

For Section 4, 5 and 6, be brief and only disclose enough information for experts to understand the data analysis. For the numerical details, refer the readers to the code base on GitHub.

Modify the correspondance between the SI and the main accordingly if the figure orders change. 