
.reverseRows = function(Matrix)
{
  ind = seq(from=dim(Matrix)[1], to=1, by=-1);
  Matrix[ind,];
  #Matrix
}

.extend = function(x, n)
{
  nRep = ceiling(n/length(x));
  rep(x, nRep)[1:n];
}

# Adapt a numeric index to a subset
# Aim: if 'index' is a numeric index of special entries of a vector,
#    create a new index that references 'subset' elements of the vector  
.restrictIndex = function(index, subset)
{
  out = match(index, subset);
  out[!is.na(out)];
}

mylabeledHeatmap = function (
  Matrix, 
  xLabels, yLabels = NULL, 
  xSymbols = NULL, ySymbols = NULL, 
  colorLabels = NULL, 
  xColorLabels = FALSE, yColorLabels = FALSE,
  checkColorsValid = TRUE,
  invertColors = FALSE, 
  setStdMargins = TRUE,
  xLabelsPosition = "bottom",
  xLabelsAngle = 45,
  xLabelsAdj = 1,
  yLabelsPosition = "left",
  xColorWidth = 2*strheight("M"),
  yColorWidth = 2*strwidth("M"),
  xColorOffset = strheight("M")/3, 
  yColorOffset = strwidth("M")/3,
  # Content of heatmap
  colors = NULL, 
  naColor = "grey",
  textMatrix = NULL, cex.text = NULL, 
  textAdj = c(0.5, 0.5),
  # labeling of rows and columns
  cex.lab = NULL, 
  cex.lab.x = cex.lab,
  cex.lab.y = cex.lab,
  cex.legend = cex.lab,
  colors.lab.x = 1,
  colors.lab.y = 1,
  font.lab.x = 1,
  font.lab.y = 1,
  bg.lab.x = NULL,
  bg.lab.y = NULL,
  x.adj.lab.y = 1,
  plotLegend = TRUE, 
  keepLegendSpace = plotLegend,
  # Separator line specification                   
  verticalSeparator.x = NULL,
  verticalSeparator.col = 1,  
  verticalSeparator.lty = 1,
  verticalSeparator.lwd = 1,
  verticalSeparator.ext = 0,
  horizontalSeparator.y = NULL,
  horizontalSeparator.col = 1,  
  horizontalSeparator.lty = 1,
  horizontalSeparator.lwd = 1,
  horizontalSeparator.ext = 0,
  # optional restrictions on which rows and columns to actually show
  showRows = NULL,
  showCols = NULL,
  # Other arguments...
  ... ) 
{
  textFnc = match.fun("text");
  if (!is.null(colorLabels)) {xColorLabels = colorLabels; yColorLabels = colorLabels; }
  if (is.null(yLabels) & (!is.null(xLabels)) & (dim(Matrix)[1]==dim(Matrix)[2])) 
    yLabels = xLabels; 
  nCols = ncol(Matrix);
  nRows = nrow(Matrix);
  if (length(xLabels)!=nCols) 
    stop("Length of 'xLabels' must equal the number of columns in 'Matrix.'");
  if (length(yLabels)!=nRows)
    stop("Length of 'yLabels' must equal the number of rows in 'Matrix.'");
  if (is.null(showRows)) showRows = c(1:nRows);
  if (is.null(showCols)) showCols = c(1:nCols);
  nShowCols = length(showCols);
  nShowRows = length(showRows);
  if (nShowCols==0) stop("'showCols' is empty.");
  if (nShowRows==0) stop("'showRows' is empty.");
  if (checkColorsValid)
  {
    xValidColors = !is.na(match(substring(xLabels, 3), colors()));
    yValidColors = !is.na(match(substring(yLabels, 3), colors()));
  } else {
    xValidColors = rep(TRUE, length(xLabels));
    yValidColors = rep(TRUE, length(yLabels));
  }
  if (sum(xValidColors)>0) xColorLabInd = xValidColors[showCols]
  if (sum(!xValidColors)>0) xTextLabInd = !xValidColors[showCols]
  if (sum(yValidColors)>0) yColorLabInd = yValidColors[showRows]
  if (sum(!yValidColors)>0) yTextLabInd = !yValidColors[showRows]
  if (setStdMargins)
  {
    if (xColorLabels & yColorLabels)
    {
      par(mar=c(2,2,3,5)+0.2);
    } else {
      par(mar = c(7,7,3,5)+0.2);
    }
  }
  xLabels.show = xLabels[showCols];
  yLabels.show = yLabels[showRows];
  
  if (!is.null(xSymbols))
  {
    if (length(xSymbols)!=nCols)
      stop("When 'xSymbols' are given, their length must equal the number of columns in 'Matrix.'");
    xSymbols.show = xSymbols[showCols];
  } else 
    xSymbols.show = NULL;
  
  if (!is.null(ySymbols))
  {
    if (length(ySymbols)!=nRows)
      stop("When 'ySymbols' are given, their length must equal the number of rows in 'Matrix.'");
    ySymbols.show = ySymbols[showRows];
  } else 
    ySymbols.show = NULL;
  
  xLabPos = charmatch(xLabelsPosition, c("bottom", "top"));
  if (is.na(xLabPos))
    stop("Argument 'xLabelsPosition' must be (a unique abbreviation of) 'bottom', 'top'");
  
  yLabPos = charmatch(yLabelsPosition, c("left", "right"));
  if (is.na(yLabPos))
    stop("Argument 'yLabelsPosition' must be (a unique abbreviation of) 'left', 'right'");
  
  if (is.null(colors)) colors = heat.colors(30);
  if (invertColors) colors = rev(colors);
  
  labPos = .heatmapWithLegend(Matrix[showRows, showCols, drop = FALSE], 
                              signed = FALSE, colors = colors, naColor = naColor, cex.legend = cex.legend, 
                              plotLegend = plotLegend,  keepLegendSpace = keepLegendSpace, ...)
  plotbox = labPos$box;
  xmin = plotbox[1]; xmax = plotbox[2]; ymin = plotbox[3]; yrange = plotbox[4]-ymin;
  ymax = plotbox[4]; xrange = xmax - xmin;
  # The positions below are for showCols/showRows-restriceted data
  xLeft = labPos$xLeft;
  xRight = labPos$xRight;
  yTop = labPos$yTop;
  yBot = labPos$yBot;
  
  xspacing = labPos$xMid[2] - labPos$xMid[1];
  yspacing = abs(labPos$yMid[2] - labPos$yMid[1]);
  
  offsetx = .extend(xColorOffset, nCols)[showCols]
  offsety = .extend(yColorOffset, nRows)[showRows]
  xColW = xColorWidth;
  yColW = yColorWidth;
  
  # Additional angle-dependent offsets for x axis labels
  textOffsetY = strheight("M") * cos(xLabelsAngle/180 * pi);
  
  if (any(xValidColors)) offsetx = offsetx + xColW;
  if (any(yValidColors)) offsety = offsety + yColW;
  
  # Create the background for column and row labels.
  extension.left = par("mai")[2] * # left margin width in inches
    par("cxy")[1] / par("cin")[1]   # character size in user corrdinates/character size in inches
  extension.right = par("mai")[4] * # right margin width in inches
    par("cxy")[1] / par("cin")[1]   # character size in user corrdinates/character size in inches
  
  extension.bottom = par("mai")[1] * 
    par("cxy")[2] / par("cin")[2]- # character size in user corrdinates/character size in inches
    offsetx   
  extension.top = par("mai")[3] * 
    par("cxy")[2] / par("cin")[2]-   # character size in user corrdinates/character size in inches
    offsetx
  
  figureBox = par("usr");
  figXrange = figureBox[2] - figureBox[1];
  figYrange = figureBox[4] - figureBox[3];
  if (!is.null(bg.lab.x))
  {
    bg.lab.x = .extend(bg.lab.x, nCols)[showCols];
    if (xLabPos==1)
    {
      y0 = ymin;
      ext = extension.bottom;
      sign = 1;
    } else {
      y0 = ymax;
      ext = extension.top;
      sign = -1;
    }
    figureDims = par("pin");
    angle = xLabelsAngle/180*pi;
    ratio = figureDims[1]/figureDims[2] * figYrange/figXrange;
    ext.x = -sign * ext * 1/tan(angle)/ratio;
    ext.y = sign * ext * sign(sin(angle))
    
    offset = (sum(xValidColors)>0) * xColW + offsetx + textOffsetY;
    
    for (cc in 1:nShowCols)
      polygon(x = c(xLeft[cc], xLeft[cc], xLeft[cc] + ext.x, xRight[cc] + ext.x, xRight[cc], xRight[cc]),
              y = c(y0, y0-sign*offset[cc], y0-sign*offset[cc] - ext.y, y0-sign*offset[cc] - ext.y, 
                    y0-sign*offset[cc], y0), 
              border = bg.lab.x[cc], col = bg.lab.x[cc], xpd = TRUE);
  }
  
  if (!is.null(bg.lab.y))
  {
    bg.lab.y = .extend(bg.lab.y, nRows)
    reverseRows = TRUE;
    if (reverseRows) bg.lab.y = rev(bg.lab.y);
    bg.lab.y = bg.lab.y[showRows];
    
    if (yLabPos==1)
    {
      xl = xmin-extension.left;
      xr = xmin;
    } else {
      xl = xmax;
      xr = xmax + extension.right;
    }
    for (r in 1:nShowRows)
      rect(xl, yBot[r], xr, yTop[r],
           col = bg.lab.y[r], border = bg.lab.y[r], xpd = TRUE);
  }
  
  colors.lab.x = .extend(colors.lab.x, nCols)[showCols];
  font.lab.x = .extend(font.lab.x, nCols)[showCols];
  # Write out labels
  if (sum(!xValidColors)>0)
  {
    xLabYPos = if(xLabPos==1) ymin - offsetx- textOffsetY else ymax + offsetx + textOffsetY;
    if (is.null(cex.lab)) cex.lab = 1;
    mapply(textFnc, x = labPos$xMid[xTextLabInd], 
           y = xLabYPos, labels = xLabels.show[xTextLabInd],
           col = colors.lab.x[xTextLabInd],
           font = font.lab.x[xTextLabInd],
           MoreArgs = list(srt = xLabelsAngle, 
                           adj = xLabelsAdj, xpd = TRUE, cex = cex.lab.x));
  }
  if (sum(xValidColors)>0)
  {
    baseY = if (xLabPos==1) ymin-offsetx else  ymax + offsetx;
    deltaY = if (xLabPos==1) xColW else -xColW;
    rect(xleft = labPos$xMid[xColorLabInd] - xspacing/2, ybottom = baseY[xColorLabInd],
         xright = labPos$xMid[xColorLabInd] + xspacing/2, ytop = baseY[xColorLabInd] + deltaY,
         density = -1,  col = substring(xLabels.show[xColorLabInd], 3), 
         border = substring(xLabels.show[xColorLabInd], 3), xpd = TRUE)
    if (!is.null(xSymbols))
      mapply(textFnc, x = labPos$xMid[xColorLabInd], 
             y = baseY[xColorLabInd] -textOffsetY - sign(deltaY)* strwidth("M")/3, 
             labels = xSymbols.show[xColorLabInd],
             col = colors.lab.x[xColorLabInd],
             font = font.lab.x[xColorLabInd],
             MoreArgs = list( adj = xLabelsAdj, 
                              xpd = TRUE, srt = xLabelsAngle, cex = cex.lab.x));
  }
  x.adj.lab.y = .extend(x.adj.lab.y, nRows)[showRows]
  if (yLabPos==1)
  {
    marginWidth = par("mai")[2] / par("pin")[1] * xrange
  } else {
    marginWidth = par("mai")[4] / par("pin")[1] * xrange
  }
  xSpaceForYLabels = marginWidth-2*strwidth("M")/3 - ifelse(yValidColors[showRows], yColW, 0);
  xPosOfYLabels.relative = xSpaceForYLabels * (1-x.adj.lab.y) + offsety
  
  colors.lab.y = .extend(colors.lab.y, nRows)[showRows];
  font.lab.y = .extend(font.lab.y, nRows)[showRows];
  
  if (sum(!yValidColors)>0)
  {
    if (is.null(cex.lab)) cex.lab = 1;
    if (yLabPos==1)
    {
      x = xmin - strwidth("M")/3 - xPosOfYLabels.relative[yTextLabInd]
      adj = x.adj.lab.y[yTextLabInd]
    } else {
      x = xmax + strwidth("M")/3 + xPosOfYLabels.relative[yTextLabInd];
      adj = 1-x.adj.lab.y[yTextLabInd];
    }
    mapply(textFnc, y = labPos$yMid[yTextLabInd], labels = yLabels.show[yTextLabInd],
           adj = lapply(adj, c, 0.5),
           x = x,
           col = colors.lab.y[yTextLabInd],
           font = font.lab.y[yTextLabInd],
           MoreArgs = list(srt = 0, xpd = TRUE, cex = cex.lab.y));
  } 
  if (sum(yValidColors)>0)
  {
    if (yLabPos==1)
    {
      xl = xmin-offsety;
      xr = xmin-offsety + yColW;
      xtext = xmin - strwidth("M")/3 - xPosOfYLabels.relative[yColorLabInd];
      adj = x.adj.lab.y[yColorLabInd]
    } else {
      xl = xmax + offsety - yColW;
      xr = xmax + offsety;
      xtext = xmin + strwidth("M")/3 + xPosOfYLabels.relative[yColorLabInd]
      adj = 1-x.adj.lab.y[yColorLabInd];
    }
    
    rect(xleft = xl[yColorLabInd], ybottom = rev(labPos$yMid[yColorLabInd]) - yspacing/2,
         xright = xr[yColorLabInd], ytop = rev(labPos$yMid[yColorLabInd]) + yspacing/2, 
         density = -1,  col = substring(rev(yLabels.show[yColorLabInd]), 3), 
         border = substring(rev(yLabels.show[yColorLabInd]), 3), xpd = TRUE)
    #for (i in yColorLabInd)
    #{
    #  lines(c(xmin- offsetx, xmin- offsetx+yColW), y = rep(labPos$yMid[i] - yspacing/2, 2), col = i, xpd = TRUE)
    #  lines(c(xmin- offsetx, xmin- offsetx+yColW), y = rep(labPos$yMid[i] + yspacing/2, 2), col = i, xpd = TRUE)
    #}
    if (!is.null(ySymbols))
      mapply(textFnc, y = labPos$yMid[yColorLabInd], labels = ySymbols.show[yColorLabInd],
             adj = lapply(adj, c, 0.5),
             x = xtext, col = colors.lab.y[yColorLabInd], 
             font = font.lab.y[yColorLabInd],
             MoreArgs = list(srt = 0, xpd = TRUE, cex = cex.lab.y));
  }
  
  # Draw separator lines, if requested
  
  showCols.ext = c(if (1 %in% showCols) 0 else NULL, showCols);
  showCols.shift = if (0 %in% showCols.ext) 1 else 0;
  
  if (length(verticalSeparator.x) > 0)
  {
    if (any(verticalSeparator.x < 0 | verticalSeparator.x > nCols))
      stop("If given. 'verticalSeparator.x' must all be between 0 and the number of columns.");
    shownVertSep = verticalSeparator.x[ verticalSeparator.x %in% showCols.ext];
    verticalSeparator.x.show = .restrictIndex(verticalSeparator.x, showCols.ext)-showCols.shift;
    rowSepShowIndex = match(shownVertSep, verticalSeparator.x)
  } else
    verticalSeparator.x.show = NULL;
  
  if (length(verticalSeparator.x.show) > 0)
  {
    nLines = length(verticalSeparator.x);
    vs.col = .extend(verticalSeparator.col, nLines)[rowSepShowIndex];
    vs.lty = .extend(verticalSeparator.lty, nLines)[rowSepShowIndex];
    vs.lwd = .extend(verticalSeparator.lwd, nLines)[rowSepShowIndex];
    vs.ext = .extend(verticalSeparator.ext, nLines)[rowSepShowIndex];
    
    x.lines = ifelse(verticalSeparator.x.show>0, labPos$xRight[verticalSeparator.x.show], labPos$xLeft[1]);
    nLines.show = length(verticalSeparator.x.show);
    for (l in 1:nLines.show)
      lines(rep(x.lines[l], 2), c(ymin, ymax), col = vs.col[l], lty = vs.lty[l], lwd = vs.lwd[l]);
    
    angle = xLabelsAngle/180*pi;
    if (angle==0) angle = pi/2;
    if (xLabelsPosition =="bottom") 
    {
      sign = 1;
      y0 = ymin;
      ext = extension.bottom;
    } else {
      sign = -1;
      y0 = ymax;
      ext = extension.top;
    }
    figureDims = par("pin");
    ratio = figureDims[1]/figureDims[2] * figYrange/figXrange;
    ext.x = -sign * ext * 1/tan(angle)/ratio;
    ext.y = sign * ext * sign(sin(angle))
    offset = (sum(xValidColors)>0) * xColW + offsetx + textOffsetY;
    for (l in 1:nLines.show)
      lines(c(x.lines[l], x.lines[l], x.lines[l] + vs.ext[l] * ext.x[l]), 
            c(y0, y0-sign*offset[l], y0-sign*offset[l] - vs.ext[l] * ext.y[l]),  
            col = vs.col[l], lty = vs.lty[l], lwd = vs.lwd[l], xpd = TRUE);
  }
  
  showRows.ext = c(if (1 %in% showRows) 0 else NULL, showRows);
  showRows.shift = if (0 %in% showRows.ext) 1 else 0;
  
  if (length(horizontalSeparator.y) >0)
  {
    if (any(horizontalSeparator.y < 0 | horizontalSeparator.y > nRows))
      stop("If given. 'horizontalSeparator.y' must all be between 0 and the number of rows.");
    shownHorizSep = horizontalSeparator.y[ horizontalSeparator.y %in% showRows.ext];
    horizontalSeparator.y.show = .restrictIndex(horizontalSeparator.y, showRows.ext)-showRows.shift;
    rowSepShowIndex = match(shownHorizSep, horizontalSeparator.y)
  } else 
    horizontalSeparator.y.show = NULL;
  
  if (length(horizontalSeparator.y.show) > 0)
  {
    reverseRows = TRUE;
    if (reverseRows) 
    {
      horizontalSeparator.y.show = nShowRows - horizontalSeparator.y.show+1;
      y.lines = ifelse( horizontalSeparator.y.show <=nShowRows, 
                        labPos$yBot[horizontalSeparator.y.show], labPos$yTop[nShowRows]);
    } else {
      y.lines = ifelse( horizontalSeparator.y.show > 0, labPos$yBot[horizontalSeparator.y.show], labPos$yTop[1]);
    }
    nLines = length(horizontalSeparator.y);
    vs.col = .extend(horizontalSeparator.col, nLines)[rowSepShowIndex];
    vs.lty = .extend(horizontalSeparator.lty, nLines)[rowSepShowIndex];
    vs.lwd = .extend(horizontalSeparator.lwd, nLines)[rowSepShowIndex];
    vs.ext = .extend(horizontalSeparator.ext, nLines)[rowSepShowIndex];
    nLines.show = length(horizontalSeparator.y.show);
    for (l in 1:nLines.show)
    {
      if (yLabPos==1)
      {
        xl = xmin-vs.ext[l]*extension.left;
        xr = xmax;
      } else {
        xl = xmin;
        xr = xmax + vs.ext[l]*extension.right;
      }
      
      lines(c(xl, xr), rep(y.lines[l], 2), 
            col = vs.col[l], lty = vs.lty[l], lwd = vs.lwd[l], xpd = TRUE);
    }
  }
  
  if (!is.null(textMatrix))
  {
    if (is.null(cex.text)) cex.text = par("cex");
    if (is.null(dim(textMatrix)))
      if (length(textMatrix)==prod(dim(Matrix))) dim(textMatrix)=dim(Matrix);
      if (!isTRUE(all.equal(dim(textMatrix), dim(Matrix))))
        stop("labeledHeatmap: textMatrix was given, but has dimensions incompatible with Matrix.");
      for (rw in 1:nShowRows)
        for (cl in 1:nShowCols)
        {
          text(labPos$xMid[cl], labPos$yMid[rw],
               as.character(textMatrix[showRows[rw],showCols[cl]]), xpd = TRUE, cex = cex.text, adj = textAdj);
        }
  }
  axis(1, labels = FALSE, tick = FALSE)
  axis(2, labels = FALSE, tick = FALSE)
  axis(3, labels = FALSE, tick = FALSE)
  axis(4, labels = FALSE, tick = FALSE)
  invisible(labPos)
}



# Replacement for the function image.plot

.autoTicks = function(min, max, maxTicks = 6 , tickPos = c(1,2,5))
{
  range = max - min;
  tick0 = range/maxTicks;
  maxTick = max(tickPos);
  # Ticks can only be multiples of tickPos
  mult = 1;
  if (tick0 < maxTick/10)
  {
    while (tick0 < maxTick/10) {tick0 = 10*tick0; mult = mult*10; }
  } else
    while (tick0 >=maxTick ) {tick0 = tick0/10; mult = mult/10; }

  ind = sum(tick0 > tickPos) + 1;
  tickStep = tickPos[ind] / mult;

  lowTick = min/tickStep;
  if (floor(lowTick)!=lowTick) lowTick = lowTick + 1;
  lowTick = floor(lowTick);

  ticks = tickStep * (lowTick:(lowTick + maxTicks+1));
  ticks = ticks[ticks <= max];
  ticks;
}

.plotStandaloneLegend = function(
  colors,
  lim,
  ## These dimensions are in inches
  tickLen = 0.09,
  tickGap = 0.04,
  minBarWidth = 0.09,
  maxBarWidth = Inf,
  mar = c(0.5, 0.2, 0.5, 0.1))
{
  par(mar = mar);
  plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "");
  box = par("usr");
  tickVal = .autoTicks(lim[1], lim[2]);
  pin = par("pin");
  xrange = box[2] - box[1];
  tickLen.usr = tickLen/pin[1] * xrange
  tickGap.usr = tickGap/pin[1] * xrange
  minBarWidth.usr = minBarWidth/pin[1] * xrange
  maxBarWidth.usr = maxBarWidth/pin[1] * xrange
  maxTickWidth = max(strwidth(tickVal));
  if (maxTickWidth + tickLen.usr + tickGap.usr > box[2]-box[1]-minBarWidth.usr)
    warning("Some tick labels will be truncated.");
  xMax = max(box[2]-maxTickWidth - tickLen.usr - tickGap.usr, box[1] + minBarWidth.usr);
  if (xMax - box[1] > maxBarWidth.usr) xMax = box[1] + maxBarWidth.usr;
  .plotColorLegend(box[1], xMax,
                   box[3], box[4],
                   colors = colors,
                   lim = lim,
                   tickLen.usr = tickLen.usr,
                   tickGap.usr = tickGap.usr);
}

.plotColorLegend = function(xmin, xmax, ymin, ymax,
                            colors,
                            tickLen.usr = 0.7* strwidth("M"),
                            tickGap.usr = 0.3 * strwidth("M"),
                            lim, cex.legend = 1)
{
  tickVal = .autoTicks(lim[1], lim[2]);
  tickY = (tickVal - lim[1]) / (lim[2] - lim[1]) * (ymax - ymin) + ymin;
  nTicks = length(tickVal);
  
  # Ticks:
  for (t in 1:nTicks)
    lines(c(xmax, xmax + tickLen.usr), c(tickY[t], tickY[t]), xpd = TRUE);
  text(rep(xmax + tickLen.usr + tickGap.usr), tickY, tickVal, adj = c(0, 0.5), cex = cex.legend,
       xpd = TRUE);
  
  # Fill with color:
  nColors = length(colors);
  ybl = (ymax-ymin)/nColors * (0:(nColors-1)) + ymin;
  ytl = (ymax-ymin)/nColors * (1:nColors) + ymin;
  rect(xleft = rep(xmin, nColors), xright = rep(xmax, nColors),
       ybottom = ybl, ytop = ytl, col = colors, border = colors, xpd = TRUE);
  
  lines(c(xmin, xmax, xmax, xmin, xmin), c(ymin, ymin, ymax, ymax, ymin), xpd = TRUE );
}


.heatmapWithLegend = function(data, signed, colors, naColor = "grey", zlim = NULL, 
                              reverseRows = TRUE,
                              plotLegend = TRUE,
                              keepLegendSpace = plotLegend,
                              cex.legend = cex.legend, 
                              legendShrink = 0.94,
                              ## The following arguments are now in inches
                              legendSpace = 0.5,   
                              legendWidth = 0.13,
                              legendGap = 0.09,
                              frame = TRUE,
                              frameTicks = FALSE, tickLen = 0.09,
                              ...)
{
  data = as.matrix(data); nCols = ncol(data); nRows = nrow(data);
  if (is.null(zlim)) 
  {
    zlim = range(data, na.rm = TRUE);
    if (signed) zlim = c(-max(abs(zlim)), max(abs(zlim)));
  }
  
  barplot(1, col = "white", border = "white", axisnames = FALSE,
          axes = FALSE, ...);
  
  pin = par("pin");
  box = par("usr");
  xminAll = box[1]; 
  xmaxAll = box[2]; 
  yminAll = box[3]; 
  ymaxAll = box[4]; 
  
  legendSpace.usr = legendSpace/pin[1] * (xmaxAll-xminAll);
  legendWidth.usr = legendWidth/pin[1] * (xmaxAll-xminAll);
  legendGap.usr = legendGap/pin[1] * (xmaxAll-xminAll);
  tickLen.usr = tickLen/pin[1] * (xmaxAll-xminAll);
  
  if (!keepLegendSpace && !plotLegend)
  {
    legendSpace.usr = 0;
    legendWidth.usr = 0;
    legendGap.usr = 0;
  }
  
  ymin = yminAll; 
  ymax = ymaxAll; 
  xmin = xminAll; 
  xmax = xmaxAll - legendSpace.usr;
  if (xmax < xmin) stop("'legendSpace is too large, not enough space for the heatmap."); 
  
  xStep = (xmax - xmin)/nCols; 
  xLeft = xmin + c(0:(nCols-1)) * xStep;
  xRight = xLeft + xStep; 
  xMid = (xLeft + xRight)/2;
  
  yStep = (ymax - ymin)/nRows; yBot  = ymin + c(0:(nRows-1)) * yStep;
  yTop  = yBot + yStep; yMid = c(yTop+ yBot)/2;
  
  if (reverseRows)
  {
    colorMat = numbers2colors(.reverseRows(data), signed, colors = colors, lim = zlim,
                              naColor = naColor)
  } else
    colorMat = numbers2colors(data, signed, colors = colors, lim = zlim, naColor = naColor)
  
  dim(colorMat) = dim(data);
  
  for (c in 1:nCols)
  {
    rect(xleft = rep(xLeft[c], nRows), xright = rep(xRight[c], nRows),
         ybottom = yBot, ytop = yTop, col = colorMat[, c], border = colorMat[, c]);
  }
  if (frame) lines( c(xmin, xmax, xmax, xmin, xmin), c(ymin, ymin, ymax, ymax, ymin) );
  
  if (plotLegend)
  {
    # Now plot the legend.
    .plotColorLegend(xmin = xmaxAll - (legendSpace.usr - legendGap.usr),
                     xmax = xmaxAll - (legendSpace.usr - legendGap.usr - legendWidth.usr),
                     ymin = yminAll + (1-legendShrink) * (ymaxAll - yminAll),
                     ymax =  ymaxAll - (1-legendShrink) * (ymaxAll - yminAll),
                     lim = zlim,
                     colors = colors,
                     tickLen.usr = tickLen.usr,
                     cex.legend = cex.legend);
    
  }
  
  list(xMid = xMid, yMid = if (reverseRows) rev(yMid) else yMid, 
       box = c(xmin, xmax, ymin, ymax), xLeft = xLeft, xRight = xRight,
       yTop = yTop, yBot = yBot);
  
}

