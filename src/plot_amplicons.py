import sys
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.patches import Ellipse, Rectangle, Arc
import matplotlib.ticker as ticker
from matplotlib import gridspec
import numpy as np
import matplotlib.patches as mpatches

# 1. change the x and y ticks for cycle
# 2. combine the pdf output for arc and cycle graph
# 3. add arg option in code to be able to adjust

class graph:
    
    # graphfile
    # shortread
    # longread
    # maxCN
    # SeqEdge
        # [chr, StartPosition, EndPosition, PredictedCopyCount, 
        # AverageCoverage, Size, NumberReadsMapped]
    # BpEdge
        # [type, chr, StartChr, StartPosition, StartStrand, EndChr, EndPosition, 
        # EndStrand, PredictedCopyCount, NumberOfReadPairs, 
        # HomologySizeIfAvailable(<0ForInsertions), Homology/InsertionSequence]
    # AxeParam
    # AxeParamList
    # xticklist
    # xticklabel
    # chrcutoff
    # inchrcutoff
    
    def __init__(self, graphfile, shortread, longread, cyclefile):
        self.graphfile = graphfile
        self.shortread = shortread
        self.longread = longread
        self.cyclefile = cyclefile
        
    def ParseGraph(self):
        lines = self.graphfile.readlines()
        maxCN = 0
        SeqEdge = {}
        for line in lines:
            if "SequenceEdge" in line:
                continue
            elif "#hybrid" in line:
                continue
            elif "BreakpointEdge" in line:
                break
            else:
                parse = line.strip().split("\t")
                parse[0] = parse[1].split(":")[0]
                parse[1] = int(parse[1].split(":")[1][:-1])
                parse[2] = int(parse[2].split(":")[1][:-1])
                parse[3] = float(parse[3])
                parse[4] = float(parse[4])
                parse[5] = int(parse[5])
                parse[6] = int(parse[6])
                if parse[3] > maxCN:
                    maxCN = parse[3]
                if parse[0] in SeqEdge:
                    SeqEdge[parse[0]].append(parse)
                else:
                    SeqEdge[parse[0]] = []
                    SeqEdge[parse[0]].append(parse)
        self.maxCN = maxCN
        self.SeqEdge = SeqEdge

        BpEdge = []
        bp = 0
        for i in range(len(lines)):
            if "BreakpointEdge" in lines[i]:
                bp = 1
                continue
            if bp == 1:
                parse = lines[i].strip().split("\t")
                if parse[0] != "discordant":
                    continue
                else:
                    start = parse[1].split("->")[0]
                    end = parse[1].split("->")[1]
                    parse.insert(1,start.split(":")[0])
                    parse.insert(2,int(start.split(":")[1][:-1]))
                    parse.insert(3,start.split(":")[1][-1:])
                    parse.insert(4,end.split(":")[0])
                    parse.insert(5,int(end.split(":")[1][:-1]))
                    parse[6] = end.split(":")[1][-1:]
                    parse[7] = float(parse[7])
                    parse[8] = int(parse[8])
                    parse[9] = int(parse[9])
                    BpEdge.append(parse)
            else:
                continue
        shortreadsum = 0
        longreadsum = 0
        countshortread = 0   
        countlongread = 0   
        for element in BpEdge:
            if element[0] == "discordant":
                if element[8] > 0:
                    shortreadsum += element[8]
                    countshortread += 1
                if element[9] > 0:
                    longreadsum += element[9]
                    countlongread += 1
        if countshortread == 0:
            countshortread += 0.1
        if countlongread == 0:
            countlongread += 0.1
        shortcopyavg = shortreadsum/countshortread
        longcopyavg = longreadsum/countlongread
        self.shortcopyavg = shortcopyavg
        self.longcopyavg = longcopyavg
        self.BpEdge = BpEdge
    
    # m is the margin between different sections. m was originally setted to 1
    def getAxe(self, m=1):
        numchr = len(self.SeqEdge)

        # actual start and end loc of each block
        AxeParam = {}
        # list of start and end loc of each block
        AxeParamList = []

        totalrange = 0
        gapsize = 0
        widthr = []

        # get intervals from SeqEdge
        for chr in self.SeqEdge:
            min = self.SeqEdge[chr][0][1]
            max = min
            start = 0
            end = 0
            AxeParam[chr] = []
            for seq in self.SeqEdge[chr]:
                if start == 0:
                    start = seq[1]
                    end = seq[2]
                    continue
                start = seq[1]
                if start != end+1:
                    max = end
                    interval = [min,max]
                    AxeParam[chr].append(interval)
                    min = max = start
                end = seq[2]
            max = self.SeqEdge[chr][-1][2]
            interval = [min,max]
            AxeParam[chr].append(interval)


        # create AxeParamList with
        # [chr, startloc, endloc, startlocnorm, endlocnorm]
        for chr in AxeParam:
            for pair in AxeParam[chr]:
                totalrange += pair[1] - pair[0]
                paramlist = [chr, pair[0], pair[1]]
                AxeParamList.append(paramlist)

        margin = totalrange*m/100
        numinterval = len(AxeParamList)
        margintotalrange = totalrange + margin*2*numinterval

        ratio = margin*100/margintotalrange
        gapnum = 1
        for chr in AxeParam:
            for pair in AxeParam[chr]:
                gap = (pair[1]-pair[0])*100/margintotalrange
                pair.append(ratio)
                ratio += gap
                pair.append(ratio)
                ratio += 2*margin*100/margintotalrange
                pair.append(gapnum)
                gapnum += 2     

        # get xtick list and labels
        xticklist = []
        xticklabel = []

        for chr in AxeParam:
            for pair in AxeParam[chr]:
                xticklabel.append(pair[0])
                xticklabel.append(pair[1])
                xticklist.append(pair[2])
                xticklist.append(pair[3])

        # get interval cut-offs
        chrcutoff = []
        inchrcutoff = []
        for chr in AxeParam:
            chrcutoff.append(AxeParam[chr][-1][3] + margin*100/margintotalrange)
            for interv in AxeParam[chr]:
                inchrcutoff.append(interv[3] + margin*100/margintotalrange)
            inchrcutoff = inchrcutoff[:-1]
        chrcutoff = chrcutoff[:-1]
        
        # label chromosome
        # get chr label and ticks
        chrlabels = []
        chrticks = []
        for chr in AxeParam:
            chrlabel = "\n" + chr
            chrlabels.append(chrlabel)
            chrinterval = 0
            n = 0
            for axe in AxeParam[chr]:
                chrinterval += axe[3] - axe[2]
                n += 1
            chrtick = AxeParam[chr][0][2] + chrinterval/2 + (n-1)*margin*100/margintotalrange
            chrticks.append(chrtick)
        
        self.numchr = numchr
        self.AxeParam = AxeParam
        self.AxeParamList = AxeParamList
        self.xticklist = xticklist
        self.xticklabel = xticklabel
        self.chrcutoff = chrcutoff
        self.inchrcutoff = inchrcutoff
        self.chrlabels = chrlabels
        self.chrticks = chrticks
        self.margin = margin
        self.margintotalrange = margintotalrange


    def coverage(self, readfile, ax, qt, read):
        maxoverallcoverage = 0
        for chromosome in self.AxeParam:
            for loc in self.AxeParam[chromosome]:
                interval = loc[1] - loc[0]
                start = loc[0]
                end = loc[1]
                print("start is: " + str(start))
                print("end is: " + str(end))

                # cover 10000bp window
                if interval > 1000000:
                    print("access no.1")
                    maxcoverage = DrawLongReadCoverage(chromosome, read, loc, start, end, 10000, ax, qt, self.margintotalrange)
                    if maxcoverage > maxoverallcoverage:
                        maxoverallcoverage = maxcoverage

                # cover 1000bp window
                elif 1000000 >= interval > 100000:
                    print("access no.2")
                    maxcoverage = DrawLongReadCoverage(chromosome, read, loc, start, end, 1000, ax, qt, self.margintotalrange)
                    if maxcoverage > maxoverallcoverage:
                        maxoverallcoverage = maxcoverage

                # cover 150 window
                else:
                    print("access no.3")
                    maxcoverage = DrawLongReadCoverage(chromosome, read, loc, start, end, 150, ax, qt, self.margintotalrange)
                    if maxcoverage > maxoverallcoverage:
                        maxoverallcoverage = maxcoverage           

        print("maxoverallcoverage is: " + str(maxoverallcoverage))
        ax.set_ylim(0, 1.5*maxoverallcoverage)
    
    def ParseCycle(self):
        # {cycle ID: chr, startloc, endloc, strand, cycle ID, iscyclic, weight}
        cyclelist = {}
        for line in self.cyclefile:
            if "#chr" in line:
                continue
            parse = line.strip().split("\t")
            parse[1] = int(parse[1])
            parse[2] = int(parse[2])
            parse[5] = bool(parse[5])
            parse[6] = float(parse[6])
            if parse[4] in cyclelist.keys():
                cyclelist[parse[4]].append(parse)
            else:
                cyclelist[parse[4]] = []
                cyclelist[parse[4]].append(parse)
        
        for cyclenums in cyclelist:
            for segment in cyclelist[cyclenums]:
                segchr = segment[0]
                segstart = segment[1]
                segend = segment[2]
                for interval in self.AxeParam[segchr]:
                    intstart = interval[0]
                    intend = interval[1]
                    if intstart <= segstart <= intend:
                        normstart = interval[2] + (segstart-intstart)*100/self.margintotalrange
                        normend = interval[2] + (segend-intstart)*100/self.margintotalrange
                        segment.insert(3, normstart)
                        segment.insert(4, normend)

        cyclesegnum = 0
        for num in cyclelist:
            for seg in cyclelist[num]:
                if seg[7] == True:
                    cyclesegnum += 1
        
        self.cyclesegnum = cyclesegnum
        self.cyclelist = cyclelist
        # # {segment_num: [chromosome, start_loc, end_loc, start_loc_norm, end_loc_norm]}
        # segments = {}
        # # [cyclenum, copycountnum, segmentlist, Path_constraints_satisfied_list]
        # cyclelist = []
        # for line in self.cyclefile:
        #     if line[:7] == "Segment":
        #         segment = line.strip().split("\t")
        #         segment = segment[1:]
        #         segment[0] = int(segment[0])
        #         segment[2] = int(segment[2])
        #         segment[3] = int(segment[3])
        #         segments[segment[0]] = segment
        # self.cyclefile.seek(0)
        # for line in self.cyclefile:
        #     if line[:6] == "Cycle=":
        #         cycleinfo = line.strip().split(";")
        #         cycleinfo[0] = int(cycleinfo[0].split("=")[1])
        #         cycleinfo[1] = float(cycleinfo[1].split("=")[1])
        #         cycleinfo[2] = cycleinfo[2].split("=")[1].split(",")
        #         cycleinfo[3] = cycleinfo[3].split("=")[1].split(",")
        #         cyclelist.append(cycleinfo)
        # for seg in segments:
        #     axeranges = self.AxeParam[segments[seg][1]]
        #     for range in axeranges:
        #         if range[0] <= segments[seg][2] <= range[1]:
        #             segstart = segments[seg][2]
        #             segend = segments[seg][3]
        #             segstartnorm = range[2] + (segstart-range[0])*100/self.margintotalrange
        #             segendnorm = range[2] + (segend-range[0])*100/self.margintotalrange
        #             segments[seg].append(segstartnorm)
        #             segments[seg].append(segendnorm)

        # cyclesegnum = 0
        # for cycle in cyclelist:
        #     if cycle[2][0][0] != 0 and cycle[2][-1][0] != 0:
        #         cyclesegnum += len(cycle[2])*2 + 4
        
        # self.cyclesegnum = cyclesegnum
        # self.segments = segments
        # self.cyclelist = cyclelist
    
    def plot(self, name, readtype, bottommargin = 0.33, height = 4, chrfontsize = 8, chry = -0.5, tickfontsize = 8, titlefontsize = 8, dpi = 200):
        
        width = np.max([10, 2*self.numchr])
        fig, ax = plt.subplots(figsize=(width, height))
        fig.subplots_adjust(bottom=bottommargin)
        ax2 = ax.twinx()

        ratio = self.margin*100/self.margintotalrange
        ymax = 0
        for chr in self.SeqEdge:
            for axeparam in self.AxeParam[chr]:
                intervalmin = axeparam[0]
                intervalmax = axeparam[1]
                for seq in self.SeqEdge[chr]:
                    if intervalmin <= seq[1] < intervalmax:
                        xmin = ratio
                        gap = (seq[2]-seq[1])*100/self.margintotalrange
                        ratio += gap
                        xmax = ratio
                        y = seq[3]
                        if y > ymax:
                            ymax = y
                        plt.hlines(y,xmin,xmax,color="black",lw=3)
                ratio += 2*self.margin*100/self.margintotalrange

        ax2.set_ylim(0, 2*ymax)
        plt.xlim(0, 100)

        plt.xticks(self.xticklist)
        ax.set_xticklabels(self.xticklabel, rotation=90, fontsize = tickfontsize)


        for xchr in self.chrcutoff:
            plt.axvline(x = xchr, linestyle = "--")

        for xinchr in self.inchrcutoff:
            plt.axvline(x = xinchr, linestyle = ":")

        ax.set_ylabel('Coverage')
        ax.set_title(name.split('/')[-1])
        ax2.set_ylabel('CN')

        # everted = "#bc175c"
        # forward = "#00a2ad"
        # reverse = "#E58E00"
        # discordant = "#491d88"
        # colorcode = {"+-":discordant, "++":forward, "-+":everted, "--":reverse}

        # Old colors
        everted = (139/256.0, 69/256.0, 19/256.0) # (brown)
        forward = 'magenta'
        reverse = (0/256.0, 139/256.0, 139/256.0)  # (teal)
        discordant = 'red'
        colorcode = {"+-":discordant, "++":forward, "-+":everted, "--":reverse, "interchromosomal":"blue"}

        for bps in self.BpEdge:
            startchr = bps[1]
            startloc = bps[2]
            startlocnorm = FindNormLoc(startchr, startloc, self.AxeParam, self.margintotalrange)
            endchr = bps[4]
            endloc = bps[5]
            endlocnorm = FindNormLoc(endchr, endloc, self.AxeParam, self.margintotalrange)
            if startchr != endchr:
                strandtype = "interchromosomal"
            else:
                strandtype = bps[3] + bps[6]
            width = endlocnorm - startlocnorm
            middle = startlocnorm + width/2
            # short -1, longread positive (only draw in long read)
            # short positive, longread is 0 (only draw in short read)
            if readtype == "shortread" and bps[8] != -1:
                copynum = float(bps[8])
                arc = Arc((middle, 0), width, 2*ymax, angle = 0, theta1 = 0, theta2 = 180, color = colorcode[strandtype], lw=min(1.5*(copynum*3/(self.shortcopyavg*3)), 3))
                ax2.add_patch(arc)
            elif readtype == "longread" and bps[9] != 0:
                copynum = float(bps[9])
                arc = Arc((middle, 0), width, 2*ymax, angle = 0, theta1 = 0, theta2 = 180, color = colorcode[strandtype], lw=min(1.5*(copynum*3/(self.longcopyavg*3)), 3))
                ax2.add_patch(arc)

            # draw the short start and end point
            promoter = self.margin*100/self.margintotalrange
            if bps[3] == "-":
                plt.hlines(0,startlocnorm,startlocnorm+promoter,color=colorcode[strandtype],lw=3)
            if bps[6] == "-":
                plt.hlines(0,endlocnorm,endlocnorm+promoter,color=colorcode[strandtype],lw=3)
            if bps[3] == "+":
                plt.hlines(0,startlocnorm,startlocnorm-promoter,color=colorcode[strandtype],lw=3)
            if bps[6] == "+":
                plt.hlines(0,endlocnorm,endlocnorm-promoter,color=colorcode[strandtype],lw=3)

        for i in range(len(self.chrlabels)):
            plt.text(self.chrticks[i]/100, chry, self.chrlabels[i], transform = ax.transAxes, fontsize = chrfontsize, horizontalalignment='center')

        # label color code
        # legendlabel = ["positive to negative", "positive to positive", "negative to positive", "negative to negative"]
        # ax2.legend(legendlabel, labelcolor = ["darkorange", "firebrick", "slateblue", "palevioletred"])
        
        if readtype == "shortread":
            read = self.shortread
            qt = 15
            # pn_patch = mpatches.Patch(color=discordant, label='discordant(+-)')
            # pp_patch = mpatches.Patch(color=forward, label='forward(++)')
            # np_patch = mpatches.Patch(color=everted, label='everted(-+)')
            # nn_patch = mpatches.Patch(color=reverse, label='reverse(--)')
            # plt.legend(handles=[pn_patch, pp_patch, np_patch, nn_patch], fontsize = tickfontsize)
        else:
            read = self.longread
            qt = 0
            # pn_patch = mpatches.Patch(color=discordant, label='discordant(+-)')
            # pp_patch = mpatches.Patch(color=forward, label='forward(++)')
            # np_patch = mpatches.Patch(color=everted, label='everted(-+)')
            # nn_patch = mpatches.Patch(color=reverse, label='reverse(--)')
            # plt.legend(handles=[pn_patch, pp_patch, np_patch, nn_patch], fontsize = tickfontsize)
        self.coverage(readtype, ax, qt, read)

        plt.show()
        fig.savefig(name + "_" + readtype + '.png', dpi = dpi)

    def plotcycle(self, name, chrfontsize = 6, chry = -0.5, tickfontsize = 5, titlefontsize = 8, dpi = 500):
        width = np.max([10, 2*self.numchr])
        # NEED TO BE ADJUSTED
        height = 12
        fig, ax = plt.subplots(figsize=(width, height))
        fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

        """
        for xchr in self.chrcutoff:
            plt.axvline(x = xchr, linestyle = "--")
        for xinchr in self.inchrcutoff:
            plt.axvline(x = xinchr, linestyle = ":")
        """

        # draw the cycles
        startline_y = -2
        neg_ylim = 0
        cycleticks = []
        cycleticklabels = []
        for cycle in self.cyclelist:
            firststartline_y = startline_y
            neg_ylim += 6
            startline_y += -4
            if self.cyclelist[cycle][0][7] == True:
                firstline_y = startline_y
                neg_ylim += 2*len(self.cyclelist[cycle])
                for i in range(len(self.cyclelist[cycle])):
                    segstrand = self.cyclelist[cycle][i][5]
                    xstart = self.cyclelist[cycle][i][3]
                    xend = self.cyclelist[cycle][i][4]
                    seglength = xend-xstart
                    rect = patches.Rectangle((xstart, startline_y), seglength, 1, facecolor='#FFFEF7', linewidth=1, edgecolor='black')
                    ax.add_patch(rect)

                    # draw the connections
                    extension = 1.5
                    if i != len(self.cyclelist[cycle]) - 1:
                        if segstrand == "+":
                            topxmin = xend
                            topxmax = topxmin + extension
                        if segstrand == "-":
                            topxmax = xstart
                            topxmin = topxmax - extension
                        ax.hlines(y = startline_y + 0.5, xmin = topxmin, xmax = topxmax, colors="blue", lw = 1)
                
                        nextsegstrand = self.cyclelist[cycle][i+1][5]
                        nextsegstart = self.cyclelist[cycle][i+1][3]
                        nextsegend = self.cyclelist[cycle][i+1][4]
                        if nextsegstrand == "+" and segstrand == "+":
                            ax.vlines(x = topxmax, ymin = startline_y - 0.5, ymax = startline_y + 0.5, colors="blue", lw = 1)
                            ax.hlines(y = startline_y - 0.5, xmin = nextsegstart - extension, xmax = topxmax, colors="blue", lw = 1)
                            ax.vlines(x = nextsegstart - extension, ymin = startline_y - 1.5, ymax = startline_y - 0.5, colors="blue", lw = 1)
                            ax.hlines(y = startline_y - 1.5, xmin = nextsegstart - extension, xmax = nextsegstart, colors="blue", lw = 1)
                            #ax.scatter((nextsegstart - extension + topxmax)/2 + 0.8 , startline_y - 0.5 - 0.05, marker=r'$\leftarrow$', color = "blue", linewidths=0.1)
                        if nextsegstrand == "-" and segstrand == "-":
                            ax.vlines(x = topxmin, ymin = startline_y - 0.5, ymax = startline_y + 0.5, colors="blue", lw = 1)
                            ax.hlines(y = startline_y - 0.5, xmin = topxmin, xmax = nextsegend + extension, colors="blue", lw = 1)
                            #ax.scatter((nextsegend + extension + topxmax)/2 - 0.8 , startline_y - 0.5 - 0.05, marker=r'$\rightarrow$', color = "blue", linewidths=0.1)
                            ax.vlines(x = nextsegend + extension, ymin = startline_y - 1.5, ymax = startline_y - 0.5, colors="blue", lw = 1)
                            ax.hlines(y = startline_y - 1.5, xmin = nextsegend, xmax = nextsegend + extension, colors="blue", lw = 1)
                        if nextsegstrand == "+" and segstrand == "-":
                            if xstart <= nextsegstart:
                                ax.vlines(x = topxmin, ymin = startline_y - 1.5, ymax = startline_y + 0.5, colors="blue", lw = 1)
                                ax.hlines(y = startline_y - 1.5, xmin = topxmin, xmax = nextsegstart, colors="blue", lw = 1)
                                #ax.scatter((nextsegstart+topxmin)/2 - 0.8 , startline_y - 1.5 - 0.05, marker=r'$\rightarrow$', color = "blue", linewidths=0.1)
                            if xstart > nextsegstart:
                                ax.hlines(y = startline_y - 1.5, xmin = nextsegstart - extension, xmax = nextsegstart, colors="blue", lw = 1)
                                ax.vlines(x = nextsegstart - extension, ymin = startline_y - 1.5, ymax = startline_y + 0.5, colors="blue", lw = 1)
                                ax.hlines(y = startline_y + 0.5, xmin = nextsegstart - extension, xmax = topxmin, colors="blue", lw = 1)
                                #ax.scatter((nextsegend+topxmin)/2 - 0.8 , startline_y + 0.5 - 0.05, marker=r'$\leftarrow$', color = "blue", linewidths=0.1)
                        if nextsegstrand == "-" and segstrand == "+":
                            if xend >= nextsegend:
                                ax.vlines(x = topxmax, ymin = startline_y - 1.5, ymax = startline_y + 0.5, colors="blue", lw = 1)
                                ax.hlines(y = startline_y - 1.5, xmin = nextsegend, xmax = topxmax, colors="blue", lw = 1)
                                #ax.scatter((nextsegend+topxmin)/2 + 0.8 , startline_y - 1.5 - 0.05, marker=r'$\leftarrow$', color = "blue", linewidths=0.1)
                            if xend < nextsegend:
                                ax.hlines(y = startline_y - 1.5, xmin = nextsegend, xmax = nextsegend + extension, colors="blue", lw = 1)
                                #ax.scatter((topxmax + nextsegend + extension)/2 - 0.8 , startline_y + 0.5 - 0.05, marker=r'$\rightarrow$', color = "blue", linewidths=0.1)
                                ax.vlines(x = nextsegend + extension, ymin = startline_y - 1.5, ymax = startline_y + 0.5, colors="blue", lw = 1)
                                ax.hlines(y = startline_y + 0.5, xmin = topxmax, xmax = nextsegend + extension, colors="blue", lw = 1)
                    startline_y += -2                
                
                firststrand = self.cyclelist[cycle][0][5]
                laststrand = self.cyclelist[cycle][-1][5]
                firststartloc = self.cyclelist[cycle][0][3]
                firstendloc = self.cyclelist[cycle][0][4]
                laststartloc = self.cyclelist[cycle][-1][3]
                lastendloc = self.cyclelist[cycle][-1][4]
                
                # STILL NEED TO SOLVE THE EDGE CASES
                
                if firststrand == "+" and laststrand == "-":
                    ax.hlines(y = startline_y + 2.5, xmin = 2, xmax = laststartloc, colors="blue", lw = 1)
                    ax.vlines(x = 2, ymin = startline_y + 2.5, ymax = firststartline_y - 3.5, colors="blue", lw = 1)
                    ax.hlines(y = firststartline_y - 3.5, xmin = 2, xmax = firststartloc, colors="blue", lw = 1)
                elif firststrand == "+" and laststrand == "+":
                    ax.hlines(y = startline_y + 2.5, xmin = lastendloc, xmax = lastendloc + 2, colors="blue", lw = 1)
                    ax.vlines(x = lastendloc + 2, ymin = startline_y + 1.5, ymax = startline_y + 2.5, colors="blue", lw = 1)
                    ax.hlines(y = startline_y + 1.5, xmin = firststartloc - 2, xmax = lastendloc + 2, colors="blue", lw = 1)
                    ax.vlines(x = firststartloc - 2, ymin = startline_y + 1.5, ymax = firststartline_y - 3.5, colors="blue", lw = 1)
                    ax.hlines(y = firststartline_y - 3.5, xmin = firststartloc - 2, xmax = firststartloc, colors="blue", lw = 1)
                elif firststrand == "-" and laststrand == "+":
                    ax.hlines(y = startline_y + 1.5, xmin = 2, xmax = lastendloc + 2, colors="blue", lw = 1)
                    ax.vlines(x = 2, ymin = startline_y + 1.5, ymax = firststartline_y - 3.5, colors="blue", lw = 1)
                    ax.hlines(y = firststartline_y - 3.5, xmin = 2, xmax = firstendloc + 2, colors="blue", lw = 1)
                    ax.vlines(x = firstendloc + 2, ymin = firststartline_y - 1.5, ymax = firststartline_y - 0.5, colors="blue", lw = 1)
                    ax.vlines(x = lastendloc + 2, ymin = startline_y + 1.5, ymax = startline_y + 0.5, colors="blue", lw = 1)
                



                startline_y += -4
                ax.hlines(y = startline_y, xmin = 0, xmax = 100, colors="black")
                startline_y += -2
                lastline_y = startline_y
                cycleticks.append((lastline_y + firstline_y)/2)
                cycleticklabels.append("cycle " + str(cycle) + ": CN=" + str(round(self.cyclelist[cycle][0][8], 2)))   
                neg_ylim += 4
                # for i in range(len(cycle[2])):
                #     # draw the bar
                #     segnum = int(cycle[2][i][:-1])
                #     segstrand = cycle[2][i][-1]
                #     xstart = self.segments[segnum][4]
                #     xend = self.segments[segnum][5]
                #     seglength = xend - xstart
                #     rect = patches.Rectangle((xstart, startline_y), seglength, 1, facecolor='#FF9EAA', linewidth=0.2, edgecolor='black')
                #     ax.add_patch(rect)
                #     print(str(segnum) + ": " + "start at " + str(xstart) + " end at " + str(xend))
                #     # draw the connections
                #     extension = 1.5
                #     if i != len(cycle[2]) - 1:
                #         if segstrand == "+":
                #             topxmin = xend
                #             topxmax = topxmin + extension
                #         if segstrand == "-":
                #             topxmax = xstart
                #             topxmin = topxmax - extension
                #         ax.hlines(y = startline_y + 0.5, xmin = topxmin, xmax = topxmax, colors="blue", lw = 0.2)

                    #     nextsegstrand = cycle[2][i + 1][-1]
                    #     nextsegnum = int(cycle[2][i + 1][:-1])
                    #     nextsegstart = self.segments[nextsegnum][4]
                    #     nextsegend = self.segments[nextsegnum][5]
                    #     if nextsegstrand == "+" and segstrand == "+":
                    #         ax.vlines(x = topxmax, ymin = startline_y - 0.5, ymax = startline_y + 0.5, colors="blue", lw = 0.2)
                    #         ax.hlines(y = startline_y - 0.5, xmin = nextsegstart - extension, xmax = topxmax, colors="blue", lw = 0.2)
                    #         ax.vlines(x = nextsegstart - extension, ymin = startline_y - 1.5, ymax = startline_y - 0.5, colors="blue", lw = 0.2)
                    #         ax.hlines(y = startline_y - 1.5, xmin = nextsegstart - extension, xmax = nextsegstart, colors="blue", lw = 0.2)
                    #     if nextsegstrand == "-" and segstrand == "-":
                    #         ax.vlines(x = topxmin, ymin = startline_y - 0.5, ymax = startline_y + 0.5, colors="blue", lw = 0.2)
                    #         ax.hlines(y = startline_y - 0.5, xmin = topxmin, xmax = nextsegend + extension, colors="blue", lw = 0.2)
                    #         ax.vlines(x = nextsegend + extension, ymin = startline_y - 1.5, ymax = startline_y - 0.5, colors="blue", lw = 0.2)
                    #         ax.hlines(y = startline_y - 1.5, xmin = nextsegend, xmax = nextsegend + extension, colors="blue", lw = 0.2)
                    #     if nextsegstrand == "+" and segstrand == "-":
                    #         if xstart <= nextsegstart:
                    #             ax.vlines(x = topxmin, ymin = startline_y - 1.5, ymax = startline_y + 0.5, colors="blue", lw = 0.2)
                    #             ax.hlines(y = startline_y - 1.5, xmin = topxmin, xmax = nextsegstart, colors="blue", lw = 0.2)
                    #         if xstart > nextsegstart:
                    #             ax.hlines(y = startline_y - 1.5, xmin = nextsegstart - extension, xmax = nextsegstart, colors="blue", lw = 0.2)
                    #             ax.vlines(x = nextsegstart - extension, ymin = startline_y - 1.5, ymax = startline_y + 0.5, colors="blue", lw = 0.2)
                    #             ax.hlines(y = startline_y + 0.5, xmin = nextsegstart - extension, xmax = topxmin, colors="blue", lw = 0.2)
                    #     if nextsegstrand == "-" and segstrand == "+":
                    #         if xend >= nextsegend:
                    #             ax.vlines(x = topxmax, ymin = startline_y - 1.5, ymax = startline_y + 0.5, colors="blue", lw = 0.2)
                    #             ax.hlines(y = startline_y - 1.5, xmin = nextsegend, xmax = topxmax, colors="blue", lw = 0.2)
                    #         if xend < nextsegend:
                    #             ax.hlines(y = startline_y - 1.5, xmin = nextsegend, xmax = nextsegend + extension, colors="blue", lw = 0.2)
                    #             ax.vlines(x = nextsegend + extension, ymin = startline_y - 1.5, ymax = startline_y + 0.5, colors="blue", lw = 0.2)
                    #             ax.hlines(y = startline_y + 0.5, xmin = topxmax, xmax = nextsegend + extension, colors="blue", lw = 0.2)
                    # startline_y += -2
                
                # firststrand = cycle[2][0][-1]
                # laststrand = cycle[2][-1][-1]
                # firstsegnum = int(cycle[2][0][:-1])
                # lastsegnum = int(cycle[2][-1][:-1])
                # firststartloc = self.segments[firstsegnum][4]
                # firstendloc = self.segments[firstsegnum][5]
                # laststartloc = self.segments[lastsegnum][4]
                # lastendloc = self.segments[lastsegnum][5]
                # # STILL NEED TO SOLVE THE EDGE CASES
                # if firststrand == "+" and laststrand == "-":
                #     ax.hlines(y = startline_y + 2.5, xmin = 2, xmax = laststartloc, colors="blue", lw = 0.2)
                #     ax.vlines(x = 2, ymin = startline_y + 2.5, ymax = firststartline_y + 0.5, colors="blue", lw = 0.2)
                #     ax.hlines(y = firststartline_y + 0.5, xmin = 2, xmax = firststartloc, colors="blue", lw = 0.2)
                # if firststrand == "+" and laststrand == "+":
                #     ax.hlines(y = startline_y + 2.5, xmin = lastendloc, xmax = lastendloc + 2, colors="blue", lw = 0.2)
                #     ax.vlines(x = lastendloc + 2, ymin = startline_y + 1.5, ymax = startline_y + 2.5, colors="blue", lw = 0.2)
                #     ax.hlines(y = startline_y + 1.5, xmin = firststartloc - 2, xmax = lastendloc + 2, colors="blue", lw = 0.2)
                #     ax.vlines(x = firststartloc - 2, ymin = startline_y + 1.5, ymax = firststartline_y + 0.5, colors="blue", lw = 0.2)
                #     ax.hlines(y = firststartline_y + 0.5, xmin = firststartloc - 2, xmax = firststartloc, colors="blue", lw = 0.2)                   

                # startline_y += -2
                # ax.hlines(y = startline_y, xmin = 0, xmax = 100, colors="black")
                # startline_y += -2
                # lastline_y = startline_y
                # cycleticks.append((lastline_y - firstline_y)/2)
                # cycleticklabels.append("cycle " + str(cycle[0]) + ": CN=" + str(round(cycle[1], 2)))
        
        ax.set_title(name + " cycle", fontsize = titlefontsize)
        plt.xlim(-10, 100)
        plt.ylim(-1*neg_ylim, 0)
        plt.yticks(cycleticks)
        ax.set_yticklabels(cycleticklabels, fontsize = 5)
        fig.savefig(name + "_cycle.png", dpi = dpi)



def FindNormLoc(chr, loc, AxeParam, margintotalrange):
    if len(AxeParam[chr])==1:
        numgap = AxeParam[chr][0][4]
        normloc = AxeParam[chr][0][2] + (loc - AxeParam[chr][0][0])*100/margintotalrange
    else:
        for param in AxeParam[chr]:
            intervalmax = param[1]
            intervalmin = param[0]
            if intervalmin <= loc <= intervalmax:
                numgap = param[4]
                normloc = param[2] + (loc - intervalmin)*100/margintotalrange
    return normloc

def DrawLongReadCoverage(chromosome, read, loclist, start, end, windowsize, ax, qt, margintotalrange):
    maxcoverage = 0
    starttick = loclist[2]
    endtick = loclist[3]
    normwindow = windowsize*100/margintotalrange
    print("normwindow is: " + str(normwindow))
    while start + windowsize <= end:
        coveragearray = read.count_coverage(chromosome, start, start + windowsize, quality_threshold = qt, read_callback = 'nofilter')
        coverage = sum(coveragearray[0]) + sum(coveragearray[1]) + sum(coveragearray[2]) + sum(coveragearray[3])
        coverage = coverage/windowsize
        if coverage > maxcoverage:
            maxcoverage = coverage
        rect = matplotlib.patches.Rectangle((starttick, 0),normwindow, coverage, color ='silver')
        ax.add_patch(rect)
        print("coverage is: " + str(coverage))
        print("starttick is: " + str(starttick))
        starttick += normwindow
        start += windowsize
    coveragearray = read.count_coverage(chromosome, start, end, quality_threshold = 0, read_callback = 'nofilter')
    coverage = sum(coveragearray[0]) + sum(coveragearray[1]) + sum(coveragearray[2]) + sum(coveragearray[3])
    coverage = coverage/(end-start)
    if coverage > maxcoverage:
        maxcoverage = coverage
    rect = matplotlib.patches.Rectangle((starttick, 0),endtick-starttick, coverage, color ='silver')
    ax.add_patch(rect)
    print("coverage is: " + str(coverage))
    print("starttick is: " + str(starttick))
    return maxcoverage

if __name__ == '__main__':
    graphname = sys.argv[1]
    ampname = graphname[:-10]
    longread = pysam.AlignmentFile(sys.argv[2], "rb")
    shortread = pysam.AlignmentFile(sys.argv[3], "rb")
    file = open(graphname, "r")
    cyclefile = open(graphname.replace("graph.txt", "cycles.bed"), "r")
    graph1 = graph(file, shortread, longread, cyclefile)
    graph1.ParseGraph()
    graph1.getAxe()
    graph1.ParseCycle()
    for num in range(int(sys.argv[4])):
        if graph1.cyclelist[str(num + 1)][0][7] == True:
            graph1.plotcycle(ampname)
            break
    file.close()
    cyclefile.close()
    
    # graph1.plot(ampname, "shortread")
    #graph1.plot(ampname, "longread")
    longread.close()
    shortread.close()
    print("finish")


