import matplotlib.pyplot as plt


def bar_plot(data):
    label = "Entire qatari dataset"
    width = 1
    diabetic_data = list(data.values()).count('diabetic')
    non_diabetic_data = list(data.values()).count('non-diabetic')

    canvas_size_x = 10
    canvas_size_y = 10
    fig, ax = plt.subplots(figsize=(canvas_size_x, canvas_size_y))

    ax.bar(label, diabetic_data, width=width, label=f'Diabetic (N={diabetic_data})')
    ax.bar(label, non_diabetic_data, width=width, bottom=diabetic_data, label=f'Non-Diabetic (N={non_diabetic_data})')

    ax.set_ylabel('Number of patients')
    ax.set_title(f'Cases of Type 2 Diabetes in the Qatari Genome Dataset (N = {diabetic_data + non_diabetic_data})')
    ax.set_xlim(-2, 2)
    ax.legend()
    plt.show()
