import sys
import re
import os
import pandas as pd

def parse_semester_arg(sem_arg):
    """
    Parse a string like 'spring2025' into a numeric rank: year * 10 + {spring=1, summer=2, fall=3}.
    For example, 'spring2025' -> 20251, 'fall2024' -> 20243.
    """
    pattern = r"(spring|summer|fall)(\d{4})"
    match = re.match(pattern, sem_arg.lower())
    if not match:
        raise ValueError(f"Invalid semester argument: {sem_arg}. Must be like 'spring2025' or 'fall2023'.")
    semester_str = match.group(1)   # 'spring', 'summer', or 'fall'
    year_str = match.group(2)       # '2025' (for example)
    year = int(year_str)

    semester_map = {"spring": 1, "summer": 2, "fall": 3}
    sem_rank = semester_map[semester_str]

    return year * 10 + sem_rank


def main():
    # Make sure we have the semester argument
    if len(sys.argv) < 2:
        print("Usage: python generate_members_page.py <current_semester>")
        print("Example: python generate_members_page.py spring2025")
        sys.exit(1)

    # Parse the current semester argument
    current_semester_arg = sys.argv[1]
    current_rank = parse_semester_arg(current_semester_arg)

    # Path to your Excel file
    excel_path = "/Users/bust0037/Documents/Master_Plan.xlsx"
    # Output HTML file
    output_html = "/Users/bust0037/Documents/Websites/thebustalab.github.io/members.html"

    # Read the spreadsheet into a DataFrame
    df = pd.read_excel(excel_path)

    # 1) Drop rows where PhotoPath is missing (NaN)
    df.dropna(subset=["PhotoPath"], inplace=True)
    # 2) Drop rows where PhotoPath is an empty string
    df = df[df["PhotoPath"].str.strip() != ""]

    # Map each semester to a numeric rank so we can figure out ordering
    semester_map = {"spring": 1, "summer": 2, "fall": 3}
    df["semester_rank"] = df.apply(
        lambda row: row["year"] * 10 + semester_map[row["semester"].lower()],
        axis=1
    )

    # Sort by semester_rank (ascending)
    df_sorted = df.sort_values("semester_rank")
    df_sorted = df_sorted[df_sorted['semester_rank'] <= current_rank]

    # Keep only the last (highest-rank) row for each Person
    latest_memberships = df_sorted.drop_duplicates(subset=["Person"], keep="last")
    # print(latest_memberships)

    # Also remove duplicates by PhotoPath so each image appears only once
    latest_memberships = latest_memberships.drop_duplicates(subset=["PhotoPath"], keep="last")

    # Mark who is "current" vs. "past" based on the user-supplied current_rank
    # print(current_rank)
    latest_memberships["is_current"] = latest_memberships["semester_rank"] == current_rank    

    # Split into current and past
    current_df = latest_memberships[latest_memberships["is_current"]]
    past_df = latest_memberships[~latest_memberships["is_current"]]

    # Define categorical ordering to prioritize graduate students first
    status_order = pd.Categorical(current_df["Status"], categories=["Grad", "UGRA"], ordered=True)
    current_df = current_df.assign(status_order=status_order).sort_values(
        ["status_order", "Person"], ascending=[True, True]
    )

    # Define categorical ordering to prioritize graduate students first
    status_order = pd.Categorical(past_df["Status"], categories=["Grad", "UGRA"], ordered=True)
    past_df = past_df.assign(status_order=status_order).sort_values(
        ["status_order", "Person"], ascending=[True, True]
    )

    # Hardcoded current member (Principal Investigator)
    hardcoded_member = """<div>
                            <div class="box">
                                <div class="image fit">
                                    <a href="" target="_blank">
                                        <img width=100% src="images/members/lbusta.jpg" />
                                    </a>
                                </div>
                                <div class="content">
                                        <header class="align-center">
                                            <p>Principal Investigator</p>
                                            <h3>Lucas Busta</h3>
                                            <a href="images/CV_for_website.pdf"><h4>CV</h4></a>
                                        </header>
                                </div>
                            </div>
                        </div>
"""

    # Build the HTML snippets, ensuring the hardcoded member is first in current members
    current_html = hardcoded_member + "\n" + "\n".join(generate_member_html(row) for _, row in current_df.iterrows())
    past_html    = "\n".join(generate_member_html(row) for _, row in past_df.iterrows())

    # Construct the full HTML (or just the relevant sections, if you prefer)
    full_html = f"""<!DOCTYPE HTML>
<html>
<head>
    <title>Members</title>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <link rel="stylesheet" href="assets/css/main.css" />
</head>
<body>

<!-- Header -->
    <header id="header">
        <div class="logo"><a href="index.html">busta lab<span></span></a></div>
        <a href="#menu">Menu</a>
    </header>

<!-- Nav -->
    <nav id="menu">
        <ul class="links">
            <li><a href="index.html">Home</a></li>
            <li><a href="research.html">Research</a></li>
            <li><a href="resources.html">Resources</a></li>
            <li><a href="latest_news.html">Latest News</a></li>
            <li><a href="members.html">Members</a></li>
            <li><a href="publications.html">Publications</a></li>
            <!-- <li><a href="citizen_science.html">Citizen Science</a></li> -->
            <!-- <li><a href="social_media.html">Social Media</a></li> -->
            <li><a href="contact.html">Contact</a></li>
        </ul>
    </nav>

<!-- One -->
    <section id="One" class="wrapper style3">
        <div class="inner">
            <header class="align-center">
                <p>lets work together</p>
                <h2>Members</h2>
            </header>
        </div>
    </section>

<!-- CURRENT MEMBERS SECTION -->
<section id="More" class="wrapper style2">
    <div class="inner">
        <div class="box">
            <div class="content">
                <header class="align-center">
                    <h2>Current Members</h2>
                </header>
            </div>
        </div>
        <div class="grid-style-four">
            {current_html}
        </div>
    </div>
</section>

<!-- PAST MEMBERS SECTION -->
<section id="More" class="wrapper style2">
    <div class="inner">
        <div class="box">
            <div class="content">
                <header class="align-center">
                    <h2>Past Members</h2>
                </header>
            </div>
        </div>
        <div class="grid-style-four">
            {past_html}
        </div>
    </div>
</section>

<!-- Footer -->
    <footer id="footer">
        <div class="container">
            <ul class="icons">
                <li><a href="https://twitter.com/PlantsRChemists" class="icon fa-twitter"><span class="label">Twitter</span></a></li>
                <li><a href="https://www.instagram.com/chemical_blooms/?hl=en" class="icon fa-instagram"><span class="label">Instagram</span></a></li>
                <li><a href="https://github.com/LucasBusta" class="icon fa-github"><span class="label">Github</span></a></li>
                <li><a href="mailto:bust0037@d.umn.edu" class="icon fa-envelope-o"><span class="label">Email</span></a></li>
            </ul>
        </div>
        <div class="copyright">
            &copy; Lucas Busta. All rights reserved.
        </div>
    </footer>

<!-- Scripts -->
    <script src="assets/js/jquery.min.js"></script>
    <script src="assets/js/jquery.scrollex.min.js"></script>
    <script src="assets/js/skel.min.js"></script>
    <script src="assets/js/util.js"></script>
    <script src="assets/js/main.js"></script>

</body>
</html>
"""

    # Write out the final HTML
    with open(output_html, "w", encoding="utf-8") as f:
        f.write(full_html)

    print(f"Successfully wrote updated members page to {output_html}")


def generate_member_html(row):
    """
    Given a row from our DataFrame, produce the HTML snippet
    for one member 'box'.
    """
    # Decide how to label them (Grad vs. UGRA)
    if row["Status"] == "Grad":
        status_text = "Graduate Researcher"
    else:
        status_text = "Undergraduate Researcher"

    # Basic fields
    name  = row["Person"]
    major = row.get("Major", "")       # or row["Major"] if guaranteed
    photo = row.get("PhotoPath", "")   # path to the image

    # If the image doesn't exist at PhotoPath, use the placeholder image
    if not os.path.exists(f"/Users/bust0037/Documents/Science/Websites/thebustalab.github.io/{photo}"):
        photo = "images/members/placeholder.jpg"

    snippet = f"""
<div>
  <div class="box">
    <div class="image fit">
      <a href="" target="_blank">
        <img width="100%" src="{photo}" alt="{name}" />
      </a>
    </div>
    <div class="content">
      <header class="align-center">
        <p>{status_text}</p>
        <h3>{name}</h3>
        <h4>{major}</h4>
      </header>
    </div>
  </div>
</div>
"""
    return snippet


if __name__ == "__main__":
    main()
