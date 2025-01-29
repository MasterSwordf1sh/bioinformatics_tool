import matplotlib.pyplot as plt
from PIL import Image

class Profile:
    def __init__(self, name, title, interests, education, expertise, contact_email):
        self.name = name
        self.title = title
        self.interests = interests
        self.education = education
        self.expertise = expertise
        self.contact_email = contact_email

    def display_profile(self):
        # Load background image
        background = Image.open("genomic_dna_background.jpg")  # Replace with your image path
        plt.figure(figsize=(10, 8))
        plt.imshow(background, aspect='auto')
        plt.axis('off')  # Hide axes

        # Add text on top of the image
        plt.text(0.5, 0.9, f"Name: {self.name}", fontsize=24, ha='center', color='white', weight='bold')
        plt.text(0.5, 0.85, f"Title: {self.title}", fontsize=20, ha='center', color='yellow')
        plt.text(0.5, 0.80, f"Interests: {', '.join(self.interests)}", fontsize=18, ha='center', color='cyan')
        plt.text(0.5, 0.75, f"Education: {', '.join(self.education)}", fontsize=18, ha='center', color='cyan')
        plt.text(0.5, 0.70, "Expertise:", fontsize=18, ha='center', color='white')

        for i, skill in enumerate(self.expertise):
            plt.text(0.5, 0.65 - i*0.05, f"- {skill}", fontsize=16, ha='center', color='lightgreen')

        plt.text(0.5, 0.65 - len(self.expertise)*0.05 - 0.05, f"Contact Email: {self.contact_email}", fontsize=16, ha='center', color='white')
        plt.text(0.5, 0.60 - len(self.expertise)*0.05 - 0.05, "Let's connect and explore the wonders of life sciences!", fontsize=16, ha='center', color='pink')

        plt.show()

if __name__ == "__main__":
    name = "Mxolisi Nene"
    title = "Multi-Omics | Machine Learning Enthusiast | Certified Data Scientist"
    interests = ["Machine Learning Enthusiast", "Microbial Genomics Expert"]
    education = ["Agricultural Scientist (BSc)", "Animal Breeding & Genetics (MSc, pending)"]
    expertise = [
        "Python programming",
        "C++ programming (g++)",
        "Linux operating system",
        "Scikit-learn for machine learning",
        "Machine learning for predictive modeling",
        "Microbial genomics for understanding microbial ecosystems",
        "Bioinformatics for analyzing biological data",
        "High-Performance Computing (HPC) on the CHPC cluster"
    ]
    contact_email = "mxolisinene4@gmail.com"

    my_profile = Profile(name, title, interests, education, expertise, contact_email)
    my_profile.display_profile()

