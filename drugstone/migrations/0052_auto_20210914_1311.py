# Generated by Django 3.0.5 on 2021-09-14 13:11

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('drugstone', '0051_auto_20210914_1251'),
    ]

    operations = [
        migrations.AlterField(
            model_name='network',
            name='id',
            field=models.CharField(max_length=28, primary_key=True, serialize=False, unique=True),
        ),
    ]
