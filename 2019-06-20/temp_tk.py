# exec(open('temp_tk.py').read())

# explore the mouse wheel with the Tkinter GUI toolkit
# Windows and Linux generate different events
# tested with Python25
import tkinter as tk
def mouse_wheel(event):
    global count
    # respond to Linux or Windows wheel event
    if event.num == 5 or event.delta == -120:
        count -= 1
    if event.num == 4 or event.delta == 120:
        count += 1
    spin.delete(0,"end")
    spin.insert(0,count)

def update_count():
    global count
    count = int(spin.get())
count = 0
root = tk.Tk()
root.title('turn mouse wheel')
root['bg'] = 'darkgreen'
spin = tk.Spinbox(root, from_=0, to=100, width=5, command=update_count)
spin.grid (column=10, row=0)


# with Windows OS
spin.bind("<MouseWheel>", mouse_wheel)
# with Linux OS
spin.bind("<Button-4>", mouse_wheel)
spin.bind("<Button-5>", mouse_wheel)
# label = tk.Label(root, font=('courier', 18, 'bold'), width=10)
# label.grid(column=1, row=0)
# label.pack(padx=40, pady=40)
root.mainloop()


# window = tk.Tk()
#
# spin = tk.Spinbox(window, from_=0, to=100, width=5)
# spin.grid (column=1, row=0)
#
# window.mainloop()
