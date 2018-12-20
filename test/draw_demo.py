import wx
import random

class MyPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, -1)
        
        self.Bind(wx.EVT_PAINT, self.on_paint)
        self.Bind(wx.EVT_SIZE, self.on_size)
        
        self.brushes = [wx.Brush(wx.Colour(random.random()*255,
                                  random.random()*255,
                                  random.random()*255)) for i in range(15)]
        self.pen = wx.TRANSPARENT_PEN
        
    def on_paint(self, event):
        #event.Skip()
        dc = wx.AutoBufferedPaintDC(self)
        
        w,h = self.GetClientSizeTuple()
        
        rectangles = [(0,i*h/15.,w,h/15.) for i in range(15)]
        
        dc.BeginDrawing()
        dc.DrawRectangleList(rectangles, self.pen, self.brushes)
        dc.EndDrawing()
        
    def on_size(self,event):
        self.Refresh()
        
app = wx.App()
frame = wx.Frame(None, -1, "")
panel = MyPanel(frame)

frame.Show()

app.MainLoop()
